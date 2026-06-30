% run_inverse_zoo_r0.m
% Fit [alpha, zoo_c_scale, r0] to EXPORTS-NA UVP observations.
%
% Extends the 2-param fit (alpha, zoo_c_scale) by adding microbial
% remineralization r0. BC is trusted as-is (bc_scale = 1).
%
% Obs: total BV at 125, 325, 475 m (100-2000 um)
% Method: fminsearch with log-space misfit + prior penalty
% Compare: 2-param result (alpha=0.10, zoo_c_sc=1, r0=0)
%          3-param result (this script)
%
% Expected runtime: ~15 min

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));
addpath(fullfile(script_dir, '..', 'data'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

cfg_base = SimulationConfig();
cfg_base.n_sections     = 30;
cfg_base.sinking_law    = 'kriest_8';
cfg_base.ds_kernel_mode = 'sinking_law';
cfg_base.r_to_rg        = 1.6;
cfg_base.alpha          = 0.10;
cfg_base.enable_coag    = true;
cfg_base.enable_disagg  = true;
cfg_base.disagg_mode    = 'operator_split';
cfg_base.disagg_dmax_A  = 9.39e-6 * 5;
cfg_base.enable_zoo     = true;
cfg_base.zoo_c          = 0.025;
cfg_base.zoo_s          = 1.3e-5;
cfg_base.zoo_p          = 0.5;
cfg_base.zoo_ic         = 7;
cfg_base.enable_microbe = false;
cfg_base.enable_mining  = true;

bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, 2:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

obs_depths = [125, 325, 475];
obs = load_uvp_obs(uvp_file, obs_depths);

% priors
prior.alpha       = 0.10;  sig_log_alpha = 1.5;   % loose: prior stickiness
prior.zoo_c_scale = 1.0;   sig_log_zoo   = 1.0;   % loose: allow zoo 3x up or down
prior.r0          = 0.01;  sig_log_r0    = 1.0;   % loose: center at 0.01/day

sigma_log = 0.5;   % obs uncertainty in natural log

cost = @(p) cost_zoo_r0(p, obs, obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base, prior, sigma_log, ...
    sig_log_alpha, sig_log_zoo, sig_log_r0);

% starting point: near prior
p0 = [0.10, 1.0, 0.01];
fprintf('\n--- Starting 3-param fminsearch [alpha, zoo_c_scale, r0] ---\n');
fprintf('p0: alpha=%.3f  zoo_c_sc=%.2f  r0=%.4f\n\n', p0(1), p0(2), p0(3));

opts = optimset('TolX',1e-3,'TolFun',0.1,'MaxFunEvals',400);
tic;
[p_fit, J_fit] = fminsearch(cost, p0, opts);
t_total = toc;

fprintf('\n--- Done in %.1f min ---\n', t_total/60);
fprintf('Best fit: alpha=%.4f  zoo_c_sc=%.3f  r0=%.5f  J=%.3f\n', ...
    p_fit(1), p_fit(2), p_fit(3), J_fit);

% forward runs for comparison
phi_prior = fwd_column_v3([0.10, 1.0, 0.0], obs_depths, col_grid, keps_day, ...
    prof, phi_bc_daily, n_days, cfg_base);
bv_prior  = sum(phi_prior, 2);

phi_fit = fwd_column_v3(p_fit, obs_depths, col_grid, keps_day, ...
    prof, phi_bc_daily, n_days, cfg_base);
bv_fit  = sum(phi_fit, 2);

fprintf('\nComparison at obs depths:\n');
fprintf('  %5s  %12s  %12s  %12s\n', 'z(m)', 'obs', 'prior', 'fit');
for id = 1:numel(obs_depths)
    fprintf('  %5d  %12.3e  %12.3e  %12.3e\n', ...
        obs_depths(id), obs.bv_total(id), bv_prior(id), bv_fit(id));
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 9],'Color','white');
hold on;
plot(obs.bv_total, obs_depths, 'ko', 'MarkerSize',5, 'MarkerFaceColor','k', ...
    'DisplayName','UVP obs');
plot(bv_prior, obs_depths, 'b--', 'LineWidth',1.2, ...
    'DisplayName', sprintf('prior (\\alpha=%.2f, zoo=1.0, r0=0)', 0.10));
plot(bv_fit, obs_depths, 'r-', 'LineWidth',1.2, ...
    'DisplayName', sprintf('fit (\\alpha=%.3f, zoo=%.2f, r0=%.4f)', ...
    p_fit(1), p_fit(2), p_fit(3)));
set(gca, 'YDir','reverse', 'XScale','log', 'FontSize',fs, 'Box','on', ...
    'YLim',[0 600], 'XLim',[5e-7 5e-5]);
xlabel('total BV 100-2000 \mum (m^3 m^{-3})', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','southeast', 'FontSize',fs, 'Box','off');
title('3-param fit: \alpha, zoo\_c, r_0', 'FontWeight','normal', 'FontSize',fs);

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
fig_path = fullfile(fig_dir, 'inverse_fit_zoo_r0.png');
exportgraphics(gcf, fig_path, 'Resolution', 200);
fprintf('\nFigure saved: %s\n', fig_path);


% ---- local cost function ----
function J = cost_zoo_r0(p, obs, obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base, prior, sigma_log, ...
    sig_alpha, sig_zoo, sig_r0)

alpha_try = p(1);  zoo_sc = p(2);  r0 = p(3);

if alpha_try<=0 || alpha_try>2 || zoo_sc<=0 || zoo_sc>20 || r0<0 || r0>1
    J = 1e6; return
end

phi_mod = fwd_column_v3(p, obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base);
bv_mod = sum(phi_mod, 2);

J_data = 0;
for id = 1:numel(obs_depths)
    if bv_mod(id) <= 0, J = 1e6; return, end
    J_data = J_data + ((log(bv_mod(id)) - log(obs.bv_total(id))) / sigma_log)^2;
end

J_prior = ((log(alpha_try)  - log(prior.alpha))       / sig_alpha)^2 + ...
          ((log(zoo_sc)      - log(prior.zoo_c_scale)) / sig_zoo)^2   + ...
          ((log(r0 + 1e-8)   - log(prior.r0))          / sig_r0)^2;

J = J_data + J_prior;
fprintf('  a=%.3f zoo=%.2f r0=%.4f  Jd=%.2f Jp=%.2f J=%.2f\n', ...
    alpha_try, zoo_sc, r0, J_data, J_prior, J);
end

% run_inverse_3param.m
% Fit [alpha, bc_scale, r0] to EXPORTS-NA UVP observations.
%
% alpha    : stickiness
% bc_scale : BC flux amplitude (1 = unchanged, <1 = reduce injection)
% r0       : microbial remineralization rate [day^-1]
%
% Obs: total BV at 125, 325, 475 m (100-2000 um)
% Method: fminsearch with log-space misfit + prior penalty
%
% Expected runtime: ~15-20 min

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
cfg_base.n_sections    = 30;
cfg_base.sinking_law   = 'kriest_8';
cfg_base.ds_kernel_mode = 'sinking_law';
cfg_base.r_to_rg       = 1.6;
cfg_base.alpha         = 0.10;
cfg_base.enable_coag   = true;
cfg_base.enable_disagg = true;
cfg_base.disagg_mode   = 'operator_split';
cfg_base.disagg_dmax_A = 9.39e-6 * 5;
cfg_base.enable_zoo    = true;
cfg_base.zoo_c         = 0.025;
cfg_base.zoo_s         = 1.3e-5;
cfg_base.zoo_p         = 0.5;
cfg_base.zoo_ic        = 7;
cfg_base.enable_microbe = false;   % v2 will set this per-call
cfg_base.enable_mining  = true;

k_plot = 2:10;
bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

obs_depths = [125, 325, 475];
obs = load_uvp_obs(uvp_file, obs_depths);

% prior
prior.alpha        = 0.10;  sigma_log_alpha = 1.5;
prior.bc_scale     = 0.5;   sigma_log_bc    = 0.7;   % center at 50% BC, loose
prior.r0           = 0.01;  sigma_log_r0    = 1.0;   % center at 0.01/day

sigma_log = 0.5;   % obs uncertainty in natural log space

cost = @(p) cost3(p, obs, obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base, prior, sigma_log, ...
    sigma_log_alpha, sigma_log_bc, sigma_log_r0);

% starting point
p0 = [0.05, 0.3, 0.02];
fprintf('\n--- Starting 3-param fminsearch ---\n');
fprintf('p0: alpha=%.3f  bc_scale=%.2f  r0=%.4f\n\n', p0(1), p0(2), p0(3));

opts = optimset('TolX',1e-3,'TolFun',0.1,'MaxFunEvals',300);
tic;
[p_fit, J_fit] = fminsearch(cost, p0, opts);
t_total = toc;

fprintf('\n--- Done in %.1f min ---\n', t_total/60);
fprintf('Best fit: alpha=%.4f  bc_scale=%.3f  r0=%.5f  J=%.3f\n', ...
    p_fit(1), p_fit(2), p_fit(3), J_fit);

% evaluate at fit and prior
phi_fit = fwd_column_v2(p_fit, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
bv_fit  = sum(phi_fit, 2);
phi_p0  = fwd_column_v2(p0, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
bv_p0   = sum(phi_p0, 2);

fprintf('\nComparison:\n  depth   obs          prior        fit\n');
for id = 1:numel(obs_depths)
    fprintf('  %3d m   %.3e   %.3e   %.3e\n', ...
        obs_depths(id), obs.bv_total(id), bv_p0(id), bv_fit(id));
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 9],'Color','white');
hold on;
plot(obs.bv_total, obs_depths, 'ko', 'MarkerSize',5, 'MarkerFaceColor','k', 'DisplayName','UVP obs');
plot(bv_p0,  obs_depths, 'b--', 'LineWidth',1.2, 'DisplayName', ...
    sprintf('start (\\alpha=%.2f, bc=%.1f, r0=%.3f)', p0(1), p0(2), p0(3)));
plot(bv_fit, obs_depths, 'r-',  'LineWidth',1.2, 'DisplayName', ...
    sprintf('fit (\\alpha=%.3f, bc=%.2f, r0=%.4f)', p_fit(1), p_fit(2), p_fit(3)));
set(gca,'YDir','reverse','FontSize',fs,'Box','on','XScale','log', ...
    'YLim',[0 600],'XLim',[5e-7 5e-5]);
xlabel('total BV 100-2000 \mum (m^3 m^{-3})', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('3-param fit: \alpha, bc\_scale, r0', 'FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'inverse_fit_3param.png'));
fprintf('\nSaved inverse_fit_3param.png\n');


% ---- local cost function ----
function J = cost3(p, obs, obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base, prior, sigma_log, ...
    sig_alpha, sig_bc, sig_r0)

alpha_try = p(1); bc_sc = p(2); r0 = p(3);

if alpha_try<=0 || alpha_try>2 || bc_sc<=0 || bc_sc>2 || r0<0 || r0>1
    J = 1e6; return
end

phi_mod = fwd_column_v2(p, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
bv_mod  = sum(phi_mod, 2);

J_data = 0;
for id = 1:numel(obs_depths)
    if bv_mod(id) <= 0, J = 1e6; return, end
    J_data = J_data + ((log(bv_mod(id)) - log(obs.bv_total(id))) / sigma_log)^2;
end

J_prior = ((log(alpha_try) - log(prior.alpha)) / sig_alpha)^2 + ...
          ((log(bc_sc)     - log(prior.bc_scale)) / sig_bc)^2 + ...
          ((log(r0+1e-8)   - log(prior.r0)) / sig_r0)^2;

J = J_data + J_prior;
fprintf('  a=%.3f bc=%.2f r0=%.4f  Jd=%.2f Jp=%.2f J=%.2f\n', ...
    alpha_try, bc_sc, r0, J_data, J_prior, J);
end

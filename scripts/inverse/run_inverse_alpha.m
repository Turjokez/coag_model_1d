% run_inverse_alpha.m
% Fit alpha and zoo_c_scale to EXPORTS-NA UVP observations.
%
% Parameters: [alpha, zoo_c_scale]
% Observations: total BV at 125, 325, 475 m (100-2000 um)
% Method: fminsearch (Nelder-Mead) with log-space misfit + prior penalty
%
% Expected runtime: ~10-15 min (80-100 forward calls x 6 s each)

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
cfg_base.enable_microbe = false;
cfg_base.enable_mining  = true;

k_plot = 2:10;
bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

obs_depths = [125, 325, 475];
obs = load_uvp_obs(uvp_file, obs_depths);

% prior
prior.alpha           = 0.10;
prior.sigma_log_alpha = 1.5;
prior.zoo_c_scale     = 1.0;
prior.sigma_log_zoo   = 1.0;

% cost function handle
cost = @(p) cost_fn_col(p, obs, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base, prior);

% starting point: [alpha, zoo_c_scale]
p0 = [0.03, 1.0];
fprintf('\n--- Starting fminsearch ---\n');
fprintf('p0: alpha=%.3f  zoo_c_scale=%.2f\n\n', p0(1), p0(2));

opts = optimset('TolX', 1e-3, 'TolFun', 0.1, 'MaxFunEvals', 200);
tic;
[p_fit, J_fit] = fminsearch(cost, p0, opts);
t_total = toc;

alpha_fit    = p_fit(1);
zoo_c_fit    = p_fit(2);

fprintf('\n--- fminsearch done in %.1f min ---\n', t_total/60);
fprintf('Best fit:  alpha = %.4f   zoo_c_scale = %.3f   J = %.3f\n', ...
    alpha_fit, zoo_c_fit, J_fit);

% run forward model at best-fit params
fprintf('\nRunning forward model at best-fit params...\n');
phi_fit = fwd_column(p_fit, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
bv_fit  = sum(phi_fit, 2);

% also run at starting params for comparison
phi_p0 = fwd_column(p0, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
bv_p0  = sum(phi_p0, 2);

fprintf('\nComparison at obs depths:\n');
fprintf('  depth   obs BV       prior BV     fit BV\n');
for id = 1:numel(obs_depths)
    fprintf('  %3d m   %.3e   %.3e   %.3e\n', ...
        obs_depths(id), obs.bv_total(id), bv_p0(id), bv_fit(id));
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 9],'Color','white');
hold on;
plot(obs.bv_total, obs_depths, 'ko', 'MarkerSize', 5, 'DisplayName', 'UVP obs');
plot(bv_p0,  obs_depths, 'b--', 'LineWidth', 1.2, 'DisplayName', sprintf('prior (\\alpha=%.2f)', p0(1)));
plot(bv_fit, obs_depths, 'r-',  'LineWidth', 1.2, 'DisplayName', sprintf('fit (\\alpha=%.3f)', alpha_fit));
set(gca,'YDir','reverse','FontSize',fs,'Box','on','XScale','log', ...
    'YLim',[0 600],'XLim',[5e-7 5e-5]);
xlabel('total BV 100-2000 \mum (m^3 m^{-3})', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('inverse fit: \alpha and zoo\_c', 'FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'inverse_fit_alpha.png'));
fprintf('\nSaved inverse_fit_alpha.png\n');

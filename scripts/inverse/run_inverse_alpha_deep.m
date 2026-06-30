% run_inverse_alpha_deep.m
% Fit only alpha to deep UVP observations.
%
% Notes:
% 1. zoo_c is fixed at the base value
% 2. only 325 m and 475 m are used
% 3. 125 m is skipped because it is too close to the BC layer

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

k_plot = 2:10;
bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

obs_depths = [325, 475];
obs = load_uvp_obs(uvp_file, obs_depths);

prior.alpha           = 0.10;
prior.sigma_log_alpha = 1.5;

cost = @(a) cost_fn_alpha_only(a, obs, obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base, prior);

alpha0 = 0.03;
fprintf('\n--- Starting deep-only alpha fit ---\n');
fprintf('alpha0: %.3f   zoo_c fixed: %.3f\n\n', alpha0, cfg_base.zoo_c);

opts = optimset('TolX', 1e-3, 'TolFun', 0.1, 'MaxFunEvals', 120);
tic;
[alpha_fit, J_fit] = fminsearch(cost, alpha0, opts);
t_total = toc;

fprintf('\n--- deep-only fit done in %.1f min ---\n', t_total/60);
fprintf('Best fit:  alpha = %.4f   zoo_c_scale = 1.000   J = %.3f\n', ...
    alpha_fit, J_fit);

fprintf('\nRunning forward model at best-fit alpha...\n');
phi_fit = fwd_column([alpha_fit, 1.0], obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base);
bv_fit  = sum(phi_fit, 2);

phi_p0 = fwd_column([alpha0, 1.0], obs_depths, col_grid, keps_day, prof, ...
    phi_bc_daily, n_days, cfg_base);
bv_p0  = sum(phi_p0, 2);

fprintf('\nComparison at obs depths:\n');
fprintf('  depth   obs BV       prior BV     fit BV\n');
for id = 1:numel(obs_depths)
    fprintf('  %3d m   %.3e   %.3e   %.3e\n', ...
        obs_depths(id), obs.bv_total(id), bv_p0(id), bv_fit(id));
end

fs = 7;
figure('Units','centimeters','Position',[2 2 9 8],'Color','white');
hold on;
plot(obs.bv_total, obs_depths, 'ko', 'MarkerSize', 5, 'DisplayName', 'UVP obs');
plot(bv_p0,  obs_depths, 'b--', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('prior (\\alpha=%.2f)', alpha0));
plot(bv_fit, obs_depths, 'r-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('fit (\\alpha=%.3f)', alpha_fit));
set(gca,'YDir','reverse','FontSize',fs,'Box','on','XScale','log', ...
    'YLim',[250 550],'XLim',[5e-7 2e-5]);
xlabel('total BV 100-2000 \mum (m^3 m^{-3})', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('deep-only inverse fit', 'FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'inverse_fit_alpha_deep.png'));
fprintf('\nSaved inverse_fit_alpha_deep.png\n');

% Short note:
% This run fits only alpha.
% Zoo is fixed.
% We use only deep depths where BC effect is smaller.

% test_fwd_column.m
% Call fwd_column once with default params and print result.
% Run this first to confirm the forward model works before fitting.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));
addpath(fullfile(script_dir, '..', 'data'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

% base config (same as Example 3)
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

% load BC
k_plot = 2:10;
bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

obs_depths = [125, 325, 475];

% run forward model with default params: alpha=0.10, zoo_c_scale=1
params = [0.10, 1.0];
fprintf('Running fwd_column with alpha=%.2f, zoo_c_scale=%.1f ...\n', params(1), params(2));
tic;
phi_out = fwd_column(params, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base);
t_elapsed = toc;
fprintf('Done in %.1f s\n', t_elapsed);

% print total BV at each depth
for id = 1:numel(obs_depths)
    bv = sum(phi_out(id,:));
    fprintf('  %d m: total BV = %.4e m^3/m^3\n', obs_depths(id), bv);
end

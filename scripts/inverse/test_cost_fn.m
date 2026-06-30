% test_cost_fn.m
% Evaluate cost function at a few alpha values to confirm J decreases
% as alpha decreases (model is currently too high -> less coag needed).

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

% prior: alpha log-normal centered at 0.10, loose width
prior.alpha          = 0.10;
prior.sigma_log_alpha = 1.5;   % wide: spans ~0.002 to 0.5
prior.zoo_c_scale    = 1.0;
prior.sigma_log_zoo  = 1.0;    % wide: spans ~0.1 to 10x

% test 3 alpha values, zoo_c_scale fixed at 1
alpha_test = [0.10, 0.03, 0.01];
fprintf('Testing cost at 3 alpha values (zoo_c_scale=1):\n');
for ia = 1:numel(alpha_test)
    params = [alpha_test(ia), 1.0];
    J = cost_fn_col(params, obs, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base, prior);
end

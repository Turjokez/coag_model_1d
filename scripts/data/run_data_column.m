% run_data_column.m
% First data-driven 1-D model run.
%
% Steps:
%   1. Load real eps, T, S from keps_for_dave.mat
%   2. Load UVP particle size spectrum
%   3. Map UVP surface phi to model Y0
%   4. Run ColumnSimulation
%   5. Compare model final profile with observed UVP mean profile
%
% Note:
%   This is a first smoke test. It uses observed UVP only as the initial
%   surface condition, not as daily surface forcing. So the model can fall
%   below the cruise-mean UVP profile after 30 days.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
repo_root  = fullfile(script_dir, '..', '..');
addpath(script_dir);
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- model grid ---
col_grid = ColumnGrid(1000, 20);

% --- real physics profile from VMP/keps file ---
keps_file = fullfile(repo_root, 'data', 'NA', 'Turbulance', 'keps_for_dave.mat');
profile = load_keps(keps_file, col_grid.z_centers);

% --- UVP observed particles ---
uvp_file = fullfile(repo_root, 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
uvp = parse_uvp(uvp_file);

% --- config ---
cfg = SimulationConfig( ...
    'n_sections',        30, ...
    't_final',           30, ...
    'delta_t',           0.2, ...
    'sinking_law',       'kriest_8', ...
    'ds_kernel_mode',    'sinking_law', ...
    'enable_coag',       true, ...
    'enable_sinking',    true, ...
    'enable_disagg',     true, ...
    'disagg_mode',       'operator_split', ...
    'disagg_dmax_cm',    1.0, ...       % validation fallback; profile.eps is used in ColumnRHS
    'proc_substeps',     20, ...
    'enable_surface_pp', false, ...     % first run: only observed surface Y0
    'enable_zoo',        true, ...
    'zoo_c',             0.025, ...
    'zoo_s',             1.3e-5, ...
    'zoo_p',             0.5, ...
    'zoo_ic',            7, ...
    'fp_alpha_cross',    0.5, ...
    'enable_mining',     true, ...
    'mining_s',          1.3e-5, ...
    'enable_microbe',    true, ...
    'microbe_r0',        0.003);

cfg.validate();

% --- surface Y0 from UVP top 5 m ---
mapped = map_uvp_to_model(uvp, cfg, col_grid);
Y0 = mapped.Y0_surface;

% --- run model ---
fprintf('=== run_data_column ===\n');
fprintf('Initial surface phi = %.3e cm^3/cm^3\n', sum(Y0(1, :)));
fprintf('Running t = %.1f days with real eps + UVP surface Y0...\n', cfg.t_final);

sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run('Y0', Y0);

Y_final  = squeeze(out.concentrations(end, :, :));
Yfp_final = squeeze(out.fecal_concentrations(end, :, :));
Y_total = Y_final + Yfp_final;

model_bv_z = sum(Y_total, 2);

% --- observed UVP profile mapped to model bins ---
obs_phi_depth = mapped.phi_depth;
obs_bv_depth = sum(obs_phi_depth, 2);
obs_bv_on_model_z = interp1(uvp.depth_m, obs_bv_depth, col_grid.z_centers, ...
    'linear', 'extrap');
obs_bv_on_model_z = max(obs_bv_on_model_z, 0);

% --- simple diagnostics ---
total0 = sum(Y0(:));
total_final = sum(Y_total(:));
fprintf('Final total phi = %.3e cm^3/cm^3\n', total_final);
fprintf('Change from initial = %.1f %%\n', 100 * (total_final - total0) / max(total0, eps));
fprintf('Observed UVP mean phi on model grid = %.3e cm^3/cm^3\n', sum(obs_bv_on_model_z));
fprintf('CFL = %.3f\n', out.cfl);

% --- Figure 1: depth profile ---
f1 = figure;
plot(obs_bv_on_model_z, col_grid.z_centers, 'k--', 'LineWidth', 1.2); hold on;
plot(model_bv_z,        col_grid.z_centers, 'b-',  'LineWidth', 1.2);
hold off;
set(gca, 'YDir', 'reverse');
xlabel('\phi');
ylabel('depth (m)');
legend({'UVP mean', 'model day 30'}, 'Location', 'best');
title('data column: depth profile');
saveas(f1, fullfile(fig_dir, 'data_column_depth.png'));

% --- Figure 2: surface spectrum ---
f2 = figure;
plot(1:cfg.n_sections, mapped.surface_phi, 'k--', 'LineWidth', 1.2); hold on;
plot(1:cfg.n_sections, Y_total(1, :),      'b-',  'LineWidth', 1.2);
hold off;
xlabel('section');
ylabel('\phi');
legend({'UVP surface', 'model surface day 30'}, 'Location', 'best');
title('data column: surface spectrum');
saveas(f2, fullfile(fig_dir, 'data_column_surface_spectrum.png'));

% --- Figure 3: eps profile ---
f3 = figure;
semilogx(profile.eps, col_grid.z_centers, 'k-o', 'MarkerSize', 4);
set(gca, 'YDir', 'reverse');
xlabel('\epsilon  [cm^2 s^{-3}]');
ylabel('depth (m)');
title('real eps profile');
saveas(f3, fullfile(fig_dir, 'data_column_eps.png'));

fprintf('\nSaved figures:\n');
fprintf('  %s\n', fullfile(fig_dir, 'data_column_depth.png'));
fprintf('  %s\n', fullfile(fig_dir, 'data_column_surface_spectrum.png'));
fprintf('  %s\n', fullfile(fig_dir, 'data_column_eps.png'));

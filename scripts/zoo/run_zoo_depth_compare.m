% run_zoo_depth_compare
% Compare old constant zoo vs depth-varying Stemmann zoo.

clear; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z        = profile.z;

base_cfg = SimulationConfig( ...
    'n_sections', 20, 't_final', 180, 'delta_t', 1, ...
    'sinking_law', 'kriest_8', 'ds_kernel_mode', 'sinking_law', ...
    'enable_coag', true, 'enable_sinking', true, 'enable_disagg', false, ...
    'proc_substeps', 20, 'enable_surface_pp', true, ...
    'surface_pp_bin', 1, 'surface_pp_mu', 0.1, ...
    'enable_zoo', true, 'zoo_c', 1e-4, ...
    'zoo_s', 1e-4, 'zoo_p', 0.3, 'zoo_ic', 1);

% Case 1: old constant zoo
profile_noZoo = DepthProfile(profile.z, profile.T_K, profile.S, ...
    profile.rho, profile.nu, profile.eps, profile.Kz);
cfg1 = base_cfg.copy();
cfg1.zoo_Zc = 100; cfg1.zoo_Zf = 50;
out1 = ColumnSimulation(cfg1, col_grid, profile_noZoo).run();
Yf1  = squeeze(out1.concentrations(end,:,:));

% Case 2: depth-varying Stemmann zoo
cfg2 = base_cfg.copy();
cfg2.zoo_Zc = 0.307; cfg2.zoo_Zf = 0.063;
out2 = ColumnSimulation(cfg2, col_grid, profile).run();
Yf2  = squeeze(out2.concentrations(end,:,:));

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

figure;
plot(sum(Yf1,2), z, 'b-', 'LineWidth', 1.2); hold on;
plot(sum(Yf2,2), z, 'r-', 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
legend({'const zoo (Zc=100)', 'Stemmann Fig1'}, 'Location', 'best');
title('zoo depth profile - t=180 d');
saveas(gcf, fullfile(fig_dir, 'zoo_depth_compare.png'));
fprintf('saved: zoo_depth_compare.png\n');

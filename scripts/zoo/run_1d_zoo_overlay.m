% run_1d_zoo_overlay
% Overlay plots for 1-D zoo production runs at t=180 days.

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z_plot   = col_grid.z_centers;

cfg_base = SimulationConfig( ...
    'n_sections', 20, ...
    't_final', 180, ...
    'delta_t', 1, ...
    'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_coag', true, ...
    'enable_sinking', true, ...
    'enable_disagg', false, ...
    'proc_substeps', 20, ...
    'enable_surface_pp', true, ...
    'surface_pp_bin', 1, ...
    'surface_pp_rate', 3e-8);

% Case 1: no grazing
cfg1 = cfg_base.copy();
cfg1.enable_zoo = false;
sim1 = ColumnSimulation(cfg1, col_grid, profile);
out1 = sim1.run();
Y1_f = squeeze(out1.concentrations(end, :, :));

% Case 2: Zc=100, Zf=50
cfg2 = cfg_base.copy();
cfg2.enable_zoo = true;
cfg2.zoo_Zc = 100;
cfg2.zoo_Zf = 50;
cfg2.zoo_p  = 0.3;
cfg2.zoo_ic = 1;
cfg2.zoo_c  = 1e-4;
cfg2.zoo_s  = 1e-4;
sim2 = ColumnSimulation(cfg2, col_grid, profile);
out2 = sim2.run();
Y2_f = squeeze(out2.concentrations(end, :, :));

% Case 3: Zc=200, Zf=100
cfg3 = cfg_base.copy();
cfg3.enable_zoo = true;
cfg3.zoo_Zc = 200;
cfg3.zoo_Zf = 100;
cfg3.zoo_p  = 0.3;
cfg3.zoo_ic = 1;
cfg3.zoo_c  = 1e-4;
cfg3.zoo_s  = 1e-4;
sim3 = ColumnSimulation(cfg3, col_grid, profile);
out3 = sim3.run();
Y3_f = squeeze(out3.concentrations(end, :, :));

fig_dir = fullfile(repo_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Figure 1: depth profile overlay
tot1 = sum(Y1_f, 2);
tot2 = sum(Y2_f, 2);
tot3 = sum(Y3_f, 2);

f1 = figure;
plot(tot1, z_plot, 'b-', 'LineWidth', 1.2); hold on;
plot(tot2, z_plot, 'r-', 'LineWidth', 1.2);
plot(tot3, z_plot, 'k-', 'LineWidth', 1.2); hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('depth profile at t=180 d');
legend({'no grazing', 'Zc=100', 'Zc=200'}, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_zoo_overlay_depth.png'));

% Figure 2: surface spectrum overlay (layer 1)
sec = 1:cfg_base.n_sections;
s1 = Y1_f(1, :);
s2 = Y2_f(1, :);
s3 = Y3_f(1, :);

f2 = figure;
semilogy(sec, s1, 'b-', 'LineWidth', 1.2); hold on;
semilogy(sec, s2, 'r-', 'LineWidth', 1.2);
semilogy(sec, s3, 'k-', 'LineWidth', 1.2); hold off;
xlabel('section');
ylabel('biovolume');
title('surface spectrum at t=180 d');
legend({'no grazing', 'Zc=100', 'Zc=200'}, 'Location', 'best');
saveas(f2, fullfile(fig_dir, 'run_1d_zoo_overlay_surface.png'));

fprintf('Saved:\n');
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_zoo_overlay_depth.png'));
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_zoo_overlay_surface.png'));

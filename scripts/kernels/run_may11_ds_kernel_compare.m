% run_may11_ds_kernel_compare
% Compare rectilinear vs curvilinear DS kernel in 1-D column.

clear; close all; clc;

addpath('src');
repo_root = pwd;
if ~exist('SimulationConfig', 'class')
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(fullfile(repo_root, 'src')));
end

% Shared setup.
cfg_base = SimulationConfig( ...
    'n_sections', 20, ...
    't_final', 60, ...
    'delta_t', 1, ...
    'sinking_law', 'kriest_8', ...
    'enable_coag', true, ...
    'enable_disagg', false, ...
    'proc_substeps', 20);

col_grid = ColumnGrid(1000, 20);
if isprop(col_grid, 'z_mid')
    z_plot = col_grid.z_mid;
else
    z_plot = col_grid.z_centers;
end
profile = DepthProfile.typical(z_plot);

% Run 1: rectilinear DS (legacy path).
cfg_rect = copy(cfg_base);
cfg_rect.ds_kernel_mode = 'legacy';
sim_rect = ColumnSimulation(cfg_rect, col_grid, profile);
out_rect = sim_rect.run();
Y_out_rect = out_rect.concentrations;
Y_rect = squeeze(Y_out_rect(end, :, :));   % n_z x n_sec

% Re-run rectilinear with more substeps to check stability.
cfg_rect2 = SimulationConfig('n_sections', 20, 't_final', 60, 'delta_t', 1, ...
    'sinking_law', 'kriest_8', 'ds_kernel_mode', 'legacy', ...
    'enable_coag', true, 'enable_disagg', false, 'proc_substeps', 100);
sim_rect2 = ColumnSimulation(cfg_rect2, col_grid, profile);
out_rect2 = sim_rect2.run();
Y_out_rect2 = out_rect2.concentrations;
bv0_r2 = sum(Y_out_rect2(1, :, :), 'all');
bvf_r2 = sum(Y_out_rect2(end, :, :), 'all');
fprintf('rect (substeps=100) bv change: %.4f%%\n', 100 * (bvf_r2 - bv0_r2) / bv0_r2);
Y_rect2 = squeeze(Y_out_rect2(end, :, :));
bv_rect2 = sum(Y_rect2(8:13, :), 'all');
fprintf('rect (substeps=100) midwater: %.4e\n', bv_rect2);

% Run 2: curvilinear DS (sinking-law path).
cfg_curv = copy(cfg_base);
cfg_curv.ds_kernel_mode = 'sinking_law';
sim_curv = ColumnSimulation(cfg_curv, col_grid, profile);
out_curv = sim_curv.run();
Y_out_curv = out_curv.concentrations;
Y_curv = squeeze(Y_out_curv(end, :, :));   % n_z x n_sec

% Midwater biovolume compare.
bv_rect = sum(Y_rect(8:13, :), 'all');
bv_curv = sum(Y_curv(8:13, :), 'all');
diff_pct = 100 * (bv_curv - bv_rect) / bv_rect;

fprintf('rect midwater: %.4e\n', bv_rect);
fprintf('curv midwater: %.4e\n', bv_curv);
fprintf('difference: %.2f%%\n', diff_pct);

% Biovolume change for each run.
bv0_rect = sum(Y_out_rect(1, :, :), 'all');
bvf_rect = sum(Y_out_rect(end, :, :), 'all');
fprintf('rect bv change: %.4f%%\n', 100 * (bvf_rect - bv0_rect) / bv0_rect);

bv0_curv = sum(Y_out_curv(1, :, :), 'all');
bvf_curv = sum(Y_out_curv(end, :, :), 'all');
fprintf('curv bv change: %.4f%%\n', 100 * (bvf_curv - bv0_curv) / bv0_curv);

% Figure: total biovolume vs depth.
tot_rect = sum(Y_rect, 2);
tot_curv = sum(Y_curv, 2);
xmax = max([tot_rect(:); tot_curv(:)]);

fig_dir = fullfile(repo_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

figure;
subplot(1,2,1);
plot(tot_rect, z_plot, 'b-');
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('rectilinear DS');
xlim([0, xmax]);

subplot(1,2,2);
plot(tot_curv, z_plot, 'b-');
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('curvilinear DS');
xlim([0, xmax]);

saveas(gcf, fullfile(fig_dir, 'may11_ds_kernel_compare.png'));

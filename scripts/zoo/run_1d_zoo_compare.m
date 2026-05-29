% run_1d_zoo_compare
% Compare 1-D column with and without zooplankton grazing.
% Uses Stemmann 2004 style: filter feeders + flux feeders, operator-split.

clear; close all; clc;

addpath(genpath('src'));

repo_root = pwd;
if ~exist('ColumnSimulation', 'class')
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(fullfile(repo_root, 'src')));
end

% --- shared setup ---
col_grid = ColumnGrid(1000, 20);
z_plot   = col_grid.z_centers;
profile  = DepthProfile.typical(z_plot);

cfg_base = SimulationConfig( ...
    'n_sections',    20, ...
    't_final',       60, ...
    'delta_t',       1, ...
    'sinking_law',   'kriest_8', ...
    'ds_kernel_mode','sinking_law', ...
    'enable_coag',   true, ...
    'enable_sinking',true, ...
    'enable_disagg', false, ...
    'proc_substeps', 20);

% --- case 1: no grazing ---
cfg1 = cfg_base.copy();
cfg1.enable_zoo = false;

sim1  = ColumnSimulation(cfg1, col_grid, profile);
out1  = sim1.run();
Y1    = squeeze(out1.concentrations);   % n_t x n_z x n_sec
Y1_t0 = squeeze(Y1(1, :, :));          % n_z x n_sec at t=0
Y1_tf = squeeze(Y1(end, :, :));        % n_z x n_sec at t=60

bv1_0 = sum(Y1_t0, 'all');
bv1_f = sum(Y1_tf, 'all');
fprintf('No grazing  — bv t=0: %.4e   t=60: %.4e   change: %.2f%%\n', ...
    bv1_0, bv1_f, 100*(bv1_f - bv1_0)/bv1_0);
if any(Y1_tf < 0, 'all')
    fprintf('  WARNING: negatives found in no-grazing case\n');
else
    fprintf('  negatives: none\n');
end

% --- case 2: with grazing ---
cfg2 = cfg_base.copy();
cfg2.enable_zoo = true;
cfg2.zoo_Zc     = 100;
cfg2.zoo_c      = 1e-4;
cfg2.zoo_Zf     = 50;
cfg2.zoo_s      = 1e-4;
cfg2.zoo_p      = 0.3;
cfg2.zoo_ic     = 1;

sim2  = ColumnSimulation(cfg2, col_grid, profile);
out2  = sim2.run();
Y2    = squeeze(out2.concentrations);
Y2_t0 = squeeze(Y2(1, :, :));
Y2_tf = squeeze(Y2(end, :, :));

bv2_0 = sum(Y2_t0, 'all');
bv2_f = sum(Y2_tf, 'all');
fprintf('With grazing — bv t=0: %.4e   t=60: %.4e   change: %.2f%%\n', ...
    bv2_0, bv2_f, 100*(bv2_f - bv2_0)/bv2_0);
if any(Y2_tf < 0, 'all')
    fprintf('  WARNING: negatives found in grazing case\n');
else
    fprintf('  negatives: none\n');
end

% --- figure ---
tot1 = sum(Y1_tf, 2);   % n_z x 1 total biovolume at t=60
tot2 = sum(Y2_tf, 2);

xmax = max([tot1(:); tot2(:)]);

fig_dir = fullfile(repo_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

figure;
subplot(1,2,1);
plot(tot1, z_plot, 'b-');
set(gca, 'YDir', 'reverse');
xlim([0, xmax]);
xlabel('biovolume');
ylabel('depth (m)');
title('no grazing');

subplot(1,2,2);
plot(tot2, z_plot, 'r-');
set(gca, 'YDir', 'reverse');
xlim([0, xmax]);
xlabel('biovolume');
ylabel('depth (m)');
title('with grazing');

saveas(gcf, fullfile(fig_dir, 'run_1d_zoo_compare.png'));
fprintf('Figure saved to output/figures/run_1d_zoo_compare.png\n');

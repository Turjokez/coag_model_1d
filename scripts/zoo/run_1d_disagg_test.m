% run_1d_disagg_test
% Test operator-split disaggregation in the 1-D column.
% Checks: (1) no crash, (2) no negatives, (3) disagg reduces large bins
% near surface but not at depth, (4) total bv lower than no-disagg case.

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
    'proc_substeps', 20, ...
    'enable_surface_pp', true, ...
    'surface_pp_bin', 1, ...
    'surface_pp_mu', 0.1, ...
    'enable_zoo', true, ...
    'zoo_Zc', 100, ...
    'zoo_c', 1e-4, ...
    'zoo_Zf', 50, ...
    'zoo_s', 1e-4, ...
    'zoo_p', 0.3, ...
    'zoo_ic', 1);

% Case 1: no disagg
cfg1 = cfg_base.copy();
cfg1.enable_disagg = false;
sim1 = ColumnSimulation(cfg1, col_grid, profile);
out1 = sim1.run();

% Case 2: operator-split disagg
cfg2 = cfg_base.copy();
cfg2.enable_disagg = true;
cfg2.disagg_mode   = 'operator_split';
sim2 = ColumnSimulation(cfg2, col_grid, profile);
out2 = sim2.run();

Yf1 = squeeze(out1.concentrations(end, :, :));  % n_z x n_sec
Yf2 = squeeze(out2.concentrations(end, :, :));

% basic checks
neg1 = sum(Yf1(:) < 0);
neg2 = sum(Yf2(:) < 0);
tot1 = sum(Yf1, 'all');
tot2 = sum(Yf2, 'all');

fprintf('\nChecks at t=180\n');
fprintf('%-20s  %-12s  %-12s\n', 'metric', 'no disagg', 'with disagg');
fprintf('%-20s  %12d  %12d\n',   'negatives',    neg1, neg2);
fprintf('%-20s  %12.4e  %12.4e\n', 'total bv',   tot1, tot2);
fprintf('%-20s  %12s  %12.2f\n',   'change (%)',  '-', ...
    100*(tot2-tot1)/max(tot1,eps));

% large-bin check: surface layer bins 18-20
fprintf('\nSurface layer (z=25m), large bins\n');
fprintf('%-6s  %-12s  %-12s\n', 'bin', 'no disagg', 'with disagg');
for b = 16:20
    fprintf('%-6d  %12.4e  %12.4e\n', b, Yf1(1,b), Yf2(1,b));
end

% depth profiles
prof1 = sum(Yf1, 2);
prof2 = sum(Yf2, 2);

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

f1 = figure;
plot(prof1, z_plot, 'b-', 'LineWidth', 1.2); hold on;
plot(prof2, z_plot, 'r-', 'LineWidth', 1.2); hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('disagg test - depth t=180');
legend({'no disagg', 'with disagg'}, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_disagg_test_depth.png'));

% time series of total bv
t    = out1.time;
tot1_t = squeeze(sum(sum(out1.concentrations, 2), 3));
tot2_t = squeeze(sum(sum(out2.concentrations, 2), 3));

f2 = figure;
plot(t, tot1_t, 'b-', 'LineWidth', 1.2); hold on;
plot(t, tot2_t, 'r-', 'LineWidth', 1.2); hold off;
xlabel('time (day)');
ylabel('total biovolume');
title('disagg test - total bv vs time');
legend({'no disagg', 'with disagg'}, 'Location', 'best');
saveas(f2, fullfile(fig_dir, 'run_1d_disagg_test_timeseries.png'));

fprintf('\nSaved figures to %s\n', fig_dir);

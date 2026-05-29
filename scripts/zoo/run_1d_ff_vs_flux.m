% run_1d_ff_vs_flux
% 1-D filter vs flux feeder test with density-dependent surface production.
% Production mode matches the 0-D slab: dY/dt = mu * Y(1,1), mu = 0.1 day^-1.

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
    'surface_pp_mu', 0.1, ...
    'enable_zoo', true, ...
    'zoo_Zc', 100, ...
    'zoo_c', 1e-4, ...
    'zoo_Zf', 50, ...
    'zoo_s', 1e-4, ...
    'zoo_p', 0.3, ...
    'zoo_ic', 1);

labels = {'no grazing', 'filter only', 'flux only', 'both'};
n_cases = numel(labels);
n_sec = cfg_base.n_sections;
n_z = col_grid.n_z;

tot180  = zeros(n_cases, 1);
prof180 = zeros(n_z, n_cases);
surf180 = zeros(n_sec, n_cases);

% Case 1: no grazing
cfg1 = cfg_base.copy();
cfg1.enable_zoo = false;
sim1 = ColumnSimulation(cfg1, col_grid, profile);
out1 = sim1.run();
Yf = squeeze(out1.concentrations(end, :, :));
tot180(1) = sum(Yf, 'all');
prof180(:,1) = sum(Yf, 2);
surf180(:,1) = Yf(1, :)';

% Case 2: filter only
cfg2 = cfg_base.copy();
cfg2.enable_zoo = true;
cfg2.zoo_Zf = 0;
sim2 = ColumnSimulation(cfg2, col_grid, profile);
out2 = sim2.run();
Yf = squeeze(out2.concentrations(end, :, :));
tot180(2) = sum(Yf, 'all');
prof180(:,2) = sum(Yf, 2);
surf180(:,2) = Yf(1, :)';

% Case 3: flux only
cfg3 = cfg_base.copy();
cfg3.enable_zoo = true;
cfg3.zoo_Zc = 0;
sim3 = ColumnSimulation(cfg3, col_grid, profile);
out3 = sim3.run();
Yf = squeeze(out3.concentrations(end, :, :));
tot180(3) = sum(Yf, 'all');
prof180(:,3) = sum(Yf, 2);
surf180(:,3) = Yf(1, :)';

% Case 4: both
cfg4 = cfg_base.copy();
cfg4.enable_zoo = true;
sim4 = ColumnSimulation(cfg4, col_grid, profile);
out4 = sim4.run();
Yf = squeeze(out4.concentrations(end, :, :));
tot180(4) = sum(Yf, 'all');
prof180(:,4) = sum(Yf, 2);
surf180(:,4) = Yf(1, :)';

fprintf('\n%-12s  %-14s  %-12s\n', 'run', 'total_bv_t180', 'change_%');
for i = 1:n_cases
    if i == 1
        d = 0;
    else
        d = 100 * (tot180(i) - tot180(1)) / max(tot180(1), eps);
    end
    fprintf('%-12s  %14.4e  %12.2f\n', labels{i}, tot180(i), d);
end

fprintf('\nSurface bin values at t=180\n');
fprintf('%-12s  %-12s  %-12s  %-12s\n', 'run', 'bin1', 'bin10', 'bin20');
for i = 1:n_cases
    fprintf('%-12s  %12.4e  %12.4e  %12.4e\n', ...
        labels{i}, surf180(1,i), surf180(10,i), surf180(20,i));
end

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Figure 1: depth profiles
f1 = figure;
plot(prof180(:,1), z_plot, 'k-', 'LineWidth', 1.2); hold on;
plot(prof180(:,2), z_plot, 'b-', 'LineWidth', 1.2);
plot(prof180(:,3), z_plot, 'r-', 'LineWidth', 1.2);
plot(prof180(:,4), z_plot, 'm-', 'LineWidth', 1.2); hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('filter vs flux - depth t=180');
legend(labels, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_ff_vs_flux_depth.png'));

% Figure 2: surface spectrum
sec = 1:n_sec;
f2 = figure;
semilogy(sec, surf180(:,1), 'k-', 'LineWidth', 1.2); hold on;
semilogy(sec, surf180(:,2), 'b-', 'LineWidth', 1.2);
semilogy(sec, surf180(:,3), 'r-', 'LineWidth', 1.2);
semilogy(sec, surf180(:,4), 'm-', 'LineWidth', 1.2); hold off;
xlabel('section');
ylabel('biovolume');
title('filter vs flux - surface spectrum t=180');
legend(labels, 'Location', 'best');
saveas(f2, fullfile(fig_dir, 'run_1d_ff_vs_flux_spectrum.png'));

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_ff_vs_flux_depth.png'));
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_ff_vs_flux_spectrum.png'));

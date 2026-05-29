% run_1d_Zc_sweep
% 1-D zooplankton density sweep with fixed Zc:Zf = 2:1.
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
    'zoo_c', 1e-4, ...
    'zoo_s', 1e-4, ...
    'zoo_p', 0.3, ...
    'zoo_ic', 1);

Zc_list = [50, 100, 200];
Zf_list = [25, 50, 100];
labels = {'no grazing', 'Zc=50', 'Zc=100', 'Zc=200'};

n_cases = numel(labels);
n_z = col_grid.n_z;
n_sec = cfg_base.n_sections;

tot180  = zeros(n_cases, 1);
prof180 = zeros(n_z, n_cases);
surf180 = zeros(n_sec, n_cases);

% baseline
cfg0 = cfg_base.copy();
cfg0.enable_zoo = false;
sim0 = ColumnSimulation(cfg0, col_grid, profile);
out0 = sim0.run();
Yf = squeeze(out0.concentrations(end, :, :));
tot180(1) = sum(Yf, 'all');
prof180(:,1) = sum(Yf, 2);
surf180(:,1) = Yf(1, :)';

% zoo cases
for i = 1:numel(Zc_list)
    cfg = cfg_base.copy();
    cfg.enable_zoo = true;
    cfg.zoo_Zc = Zc_list(i);
    cfg.zoo_Zf = Zf_list(i);
    sim = ColumnSimulation(cfg, col_grid, profile);
    out = sim.run();
    Yf = squeeze(out.concentrations(end, :, :));
    tot180(i+1) = sum(Yf, 'all');
    prof180(:,i+1) = sum(Yf, 2);
    surf180(:,i+1) = Yf(1, :)';
end

fprintf('\n%-12s  %-8s  %-14s  %-12s\n', 'case', 'Zc', 'total_bv_t180', 'change_%');
fprintf('%-12s  %-8s  %14.4e  %12.2f\n', labels{1}, '-', tot180(1), 0);
for i = 1:numel(Zc_list)
    d = 100 * (tot180(i+1) - tot180(1)) / max(tot180(1), eps);
    fprintf('%-12s  %-8d  %14.4e  %12.2f\n', labels{i+1}, Zc_list(i), tot180(i+1), d);
end

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Figure 1: depth profiles
f1 = figure;
plot(prof180(:,1), z_plot, 'b-', 'LineWidth', 1.2); hold on;
plot(prof180(:,2), z_plot, 'g-', 'LineWidth', 1.2);
plot(prof180(:,3), z_plot, 'r-', 'LineWidth', 1.2);
plot(prof180(:,4), z_plot, 'k-', 'LineWidth', 1.2); hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('Zc sweep - depth t=180');
legend(labels, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_Zc_sweep_depth.png'));

% Figure 2: surface spectrum
sec = 1:n_sec;
f2 = figure;
semilogy(sec, surf180(:,1), 'b-', 'LineWidth', 1.2); hold on;
semilogy(sec, surf180(:,2), 'g-', 'LineWidth', 1.2);
semilogy(sec, surf180(:,3), 'r-', 'LineWidth', 1.2);
semilogy(sec, surf180(:,4), 'k-', 'LineWidth', 1.2); hold off;
xlabel('section');
ylabel('biovolume');
title('Zc sweep - surface spectrum t=180');
legend(labels, 'Location', 'best');
saveas(f2, fullfile(fig_dir, 'run_1d_Zc_sweep_spectrum.png'));

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_Zc_sweep_depth.png'));
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_Zc_sweep_spectrum.png'));

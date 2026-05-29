% run_1d_c_sweep
% 1-D clearance-rate sweep with identifiability check.
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
    'zoo_Zf', 50, ...
    'zoo_s', 1e-4, ...
    'zoo_p', 0.3, ...
    'zoo_ic', 1);

c_vals  = [0.5e-4, 1e-4, 2e-4, 4e-4];
labels  = {'no grazing', 'c/2', 'c', '2c', '4c'};
colors  = {'k-', 'b-', 'r-', 'm-', 'g-'};
n_cases = numel(labels);
n_z = col_grid.n_z;
n_sec = cfg_base.n_sections;

tot180  = zeros(n_cases, 1);
prof180 = zeros(n_z, n_cases);
surf180 = zeros(n_sec, n_cases);
cz_prod = zeros(n_cases, 1);

% baseline
cfg0 = cfg_base.copy();
cfg0.enable_zoo = false;
sim0 = ColumnSimulation(cfg0, col_grid, profile);
out0 = sim0.run();
Yf = squeeze(out0.concentrations(end, :, :));
tot180(1) = sum(Yf, 'all');
prof180(:,1) = sum(Yf, 2);
surf180(:,1) = Yf(1, :)';
cz_prod(1) = 0;

% c cases
for i = 1:numel(c_vals)
    cfg = cfg_base.copy();
    cfg.enable_zoo = true;
    cfg.zoo_c = c_vals(i);
    sim = ColumnSimulation(cfg, col_grid, profile);
    out = sim.run();
    Yf = squeeze(out.concentrations(end, :, :));
    tot180(i+1) = sum(Yf, 'all');
    prof180(:,i+1) = sum(Yf, 2);
    surf180(:,i+1) = Yf(1, :)';
    cz_prod(i+1) = cfg.zoo_c * cfg.zoo_Zc;
end

fprintf('\n%-10s  %-14s  %-14s  %-12s\n', 'case', 'c*Zc(day^-1)', 'total_bv_t180', 'change_%');
fprintf('%-10s  %14.4e  %14.4e  %12.2f\n', labels{1}, cz_prod(1), tot180(1), 0);
for i = 2:n_cases
    d = 100 * (tot180(i) - tot180(1)) / max(tot180(1), eps);
    fprintf('%-10s  %14.4e  %14.4e  %12.2f\n', labels{i}, cz_prod(i), tot180(i), d);
end

% Identifiability check:
% Case A: c=1e-4, Zc=100, Zf=50 (same as "c" case above)
% Case B: c=2e-4, Zc=50,  Zf=25 (same c*Zc and same Zc:Zf ratio)
cfgA = cfg_base.copy();
cfgA.enable_zoo = true;
cfgA.zoo_c = 1e-4;
cfgA.zoo_Zc = 100;
cfgA.zoo_Zf = 50;
simA = ColumnSimulation(cfgA, col_grid, profile);
outA = simA.run();
totA = sum(squeeze(outA.concentrations(end, :, :)), 'all');

cfgB = cfg_base.copy();
cfgB.enable_zoo = true;
cfgB.zoo_c = 2e-4;
cfgB.zoo_Zc = 50;
cfgB.zoo_Zf = 25;
simB = ColumnSimulation(cfgB, col_grid, profile);
outB = simB.run();
totB = sum(squeeze(outB.concentrations(end, :, :)), 'all');

fprintf('\nIdentifiability check\n');
fprintf('A: c=1e-4, Zc=100, Zf=50 -> total %.4e\n', totA);
fprintf('B: c=2e-4, Zc=50,  Zf=25 -> total %.4e\n', totB);
fprintf('Difference A vs B: %.3e\n', abs(totA - totB));

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Figure 1: depth profiles
f1 = figure;
for i = 1:n_cases
    plot(prof180(:,i), z_plot, colors{i}, 'LineWidth', 1.2); hold on;
end
hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('c sweep - depth t=180');
legend(labels, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_c_sweep_depth.png'));

% Figure 2: surface spectrum
sec = 1:n_sec;
f2 = figure;
for i = 1:n_cases
    semilogy(sec, surf180(:,i), colors{i}, 'LineWidth', 1.2); hold on;
end
hold off;
xlabel('section');
ylabel('biovolume');
title('c sweep - surface spectrum t=180');
legend(labels, 'Location', 'best');
saveas(f2, fullfile(fig_dir, 'run_1d_c_sweep_spectrum.png'));

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_c_sweep_depth.png'));
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_c_sweep_spectrum.png'));

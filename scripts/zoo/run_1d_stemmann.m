% run_1d_stemmann
% Stemmann 2004 Table 3 validation in the 1-D column model.
% Same 5 cases as in the 0-D report (report_may14_zooplankton.md, Section 11).
% Checks that ZooplanktonGrazing behaves identically in the 1-D column
% when depth structure has not yet developed (10-day run, no production).
%
% Cases: no grazing, flux only, filter only, both, double Z.

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z_plot   = col_grid.z_centers;

% Stemmann 2004 Table 3 parameters.
Zc_stem = 0.307;       % ind m^-3
c_stem  = 0.025;       % m^3 ind^-1 day^-1
Zf_stem = 250;         % ind m^-3
s_stem  = 1.3e-5;      % m^2 ind^-1
p_stem  = 0.5;
ic_stem = 1;

cfg_base = SimulationConfig( ...
    'n_sections', 20, ...
    't_final', 10, ...
    'delta_t', 1, ...
    'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_coag', true, ...
    'enable_sinking', true, ...
    'enable_disagg', false, ...
    'enable_surface_pp', false, ...
    'proc_substeps', 20, ...
    'zoo_c',  c_stem, ...
    'zoo_s',  s_stem, ...
    'zoo_p',  p_stem, ...
    'zoo_ic', ic_stem);

% Five cases: no graze, flux only, filter only, both, double Z.
Zc_list = [0,       0,       Zc_stem, Zc_stem,   2*Zc_stem];
Zf_list = [0,       Zf_stem, 0,       Zf_stem,   2*Zf_stem];
labels  = {'no grazing', 'flux only', 'filter only', 'both', 'double Z'};
n_cases = numel(labels);
n_z     = col_grid.n_z;

tot_ts  = zeros(11, n_cases);   % t=0..10 days
prof_t10 = zeros(n_z, n_cases);

for i = 1:n_cases
    cfg = cfg_base.copy();
    if Zc_list(i) == 0 && Zf_list(i) == 0
        cfg.enable_zoo = false;
    else
        cfg.enable_zoo = true;
        cfg.zoo_Zc     = Zc_list(i);
        cfg.zoo_Zf     = Zf_list(i);
    end
    sim  = ColumnSimulation(cfg, col_grid, profile);
    out  = sim.run();
    Yhist = out.concentrations;   % n_t x n_z x n_sec
    for ti = 1:size(Yhist, 1)
        tot_ts(ti, i) = sum(Yhist(ti, :, :), 'all');
    end
    prof_t10(:, i) = sum(squeeze(Yhist(end, :, :)), 2);
end

t = out.time;

fprintf('\nTotal biovolume at t=10 (Stemmann 2004 Table 3 parameters)\n');
fprintf('%-14s  %-14s  %-14s\n', 'case', 'total_bv_t10', 'change_%');
bv0 = tot_ts(end, 1);
for i = 1:n_cases
    d = 100 * (tot_ts(end, i) - bv0) / max(bv0, eps);
    fprintf('%-14s  %14.4e  %12.3f\n', labels{i}, tot_ts(end, i), d);
end

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

colors = {'k-', 'r-', 'b-', 'm-', 'g--'};

% Figure 1: total biovolume vs time.
f1 = figure;
for i = 1:n_cases
    plot(t, tot_ts(:, i), colors{i}, 'LineWidth', 1.2); hold on;
end
hold off;
xlabel('time (day)');
ylabel('total biovolume');
title('Stemmann validation - total bv');
legend(labels, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_stemmann_totalvol.png'));

% Figure 2: depth profile at t=10.
f2 = figure;
for i = 1:n_cases
    plot(prof_t10(:, i), z_plot, colors{i}, 'LineWidth', 1.2); hold on;
end
hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
title('Stemmann validation - depth profile t=10');
legend(labels, 'Location', 'best');
saveas(f2, fullfile(fig_dir, 'run_1d_stemmann_depth.png'));

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_stemmann_totalvol.png'));
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_stemmann_depth.png'));

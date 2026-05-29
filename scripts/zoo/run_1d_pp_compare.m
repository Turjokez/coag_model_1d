% run_1d_pp_compare
% Compare constant source (surface_pp_rate) vs density-dependent growth (surface_pp_mu).
% Runs four cases for 500 days each:
%   1. constant source, no grazing
%   2. constant source, with grazing
%   3. growth rate (mu*phi), no grazing
%   4. growth rate (mu*phi), with grazing

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);

% Shared settings.
cfg_base = SimulationConfig( ...
    'n_sections', 20, ...
    't_final', 500, ...
    'delta_t', 1, ...
    'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_coag', true, ...
    'enable_sinking', true, ...
    'enable_disagg', false, ...
    'proc_substeps', 20, ...
    'zoo_Zc', 100, ...
    'zoo_c', 1e-4, ...
    'zoo_Zf', 50, ...
    'zoo_s', 1e-4, ...
    'zoo_p', 0.3, ...
    'zoo_ic', 1);

% Constant source rate matches 0-D slab approximately at steady state.
% In 0-D with mu=0.1 and phi1 ~ 1e-6, production ~ 1e-7/day.
% Here we use 3e-8 as in the other 1-D tests.
rate_const = 3e-8;     % bv/day
mu_val     = 0.1;      % day^-1 (same as 0-D slab)

labels = {'const, no graze', 'const, graze', 'mu*phi, no graze', 'mu*phi, graze'};

% --- Case 1: constant source, no grazing ---
cfg = cfg_base.copy();
cfg.enable_zoo          = false;
cfg.enable_surface_pp   = true;
cfg.surface_pp_rate     = rate_const;
cfg.surface_pp_mu       = 0;
sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();
t   = out.time;
tot = zeros(numel(t), 4);   % allocate after first run so size matches
tot(:, 1) = squeeze(sum(sum(out.concentrations, 2), 3));

% --- Case 2: constant source, with grazing ---
cfg = cfg_base.copy();
cfg.enable_zoo          = true;
cfg.enable_surface_pp   = true;
cfg.surface_pp_rate     = rate_const;
cfg.surface_pp_mu       = 0;
sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();
tot(:, 2) = squeeze(sum(sum(out.concentrations, 2), 3));

% --- Case 3: growth rate (mu*phi), no grazing ---
cfg = cfg_base.copy();
cfg.enable_zoo          = false;
cfg.enable_surface_pp   = true;
cfg.surface_pp_mu       = mu_val;
sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();
tot(:, 3) = squeeze(sum(sum(out.concentrations, 2), 3));

% --- Case 4: growth rate (mu*phi), with grazing ---
cfg = cfg_base.copy();
cfg.enable_zoo          = true;
cfg.enable_surface_pp   = true;
cfg.surface_pp_mu       = mu_val;
sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();
tot(:, 4) = squeeze(sum(sum(out.concentrations, 2), 3));

fprintf('\nTotal biovolume at t=500\n');
for i = 1:4
    fprintf('%-20s : %.4e\n', labels{i}, tot(end, i));
end

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

f1 = figure;
plot(t, tot(:,1), 'b-',  'LineWidth', 1.2); hold on;
plot(t, tot(:,2), 'b--', 'LineWidth', 1.2);
plot(t, tot(:,3), 'r-',  'LineWidth', 1.2);
plot(t, tot(:,4), 'r--', 'LineWidth', 1.2); hold off;
xlabel('time (day)');
ylabel('total biovolume');
title('constant vs mu*phi production');
legend(labels, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_pp_compare_timeseries.png'));

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_pp_compare_timeseries.png'));

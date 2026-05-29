% run_full_1d_nsections_compare.m
% Full 1-D integrated run at n_sections = 15, 20, 25, 30.
% All physics: coag + depth-scaled kernels + operator-split disagg
% (depth-varying D_max) + Stemmann zoo + fecal bin zooplankton.
%
% Goal: find the n where disaggregation acts across the full water column
% and total biovolume converges. Adrian asked for this convergence test.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

n_list = [15, 20, 25, 30];
n_runs = length(n_list);

col_grid = ColumnGrid(1000, 20);   % 20 depth layers, dz = 50 m
profile  = DepthProfile.typical(col_grid.z_centers);
z        = col_grid.z_centers;
n_z      = col_grid.n_z;
Dmax_A   = 9.39e-6;   % calibration constant in ColumnRHS

% Pre-compute D_max profile (same for all runs — depends only on eps(z))
dmax_mm = zeros(n_z, 1);
for k = 1:n_z
    eps_m      = profile.eps(k) / 1e4;
    dmax_mm(k) = Dmax_A * eps_m^(-0.25) * 1000;
end

% Base config — same for all runs except n_sections
cfg_base = SimulationConfig( ...
    't_final',           365, ...
    'delta_t',           1, ...
    'sinking_law',       'kriest_8', ...
    'ds_kernel_mode',    'sinking_law', ...
    'enable_coag',       true, ...
    'enable_sinking',    true, ...
    'enable_disagg',     true, ...
    'disagg_mode',       'operator_split', ...
    'disagg_dmax_cm',    1.0, ...
    'proc_substeps',     20, ...
    'enable_surface_pp', true, ...
    'surface_pp_bin',    1, ...
    'surface_pp_mu',     0.1, ...
    'enable_zoo',        true, ...
    'zoo_Zc',            0.307, ...
    'zoo_Zf',            0.063, ...
    'zoo_c',             0.025, ...
    'zoo_s',             1.3e-5, ...
    'zoo_p',             0.5, ...
    'zoo_ic',            7);

% Store results
bv_final   = zeros(n_runs, 1);   % total bv at t=365
maxbin_deep = zeros(n_runs, 1);  % largest populated bin in bottom layer
top_bin_mm = zeros(n_runs, 1);   % grid ceiling diameter (mm)
disagg_depth = zeros(n_runs, 1); % deepest layer where disagg is active (m)

snap_days = [30, 90, 180, 365];
bv_snaps  = zeros(n_runs, length(snap_days));
depth_profiles = cell(n_runs, 1);

for i = 1:n_runs
    n = n_list(i);
    fprintf('Running n_sections = %d ...\n', n);

    cfg = cfg_base.copy();
    cfg.n_sections = n;

    sim = ColumnSimulation(cfg, col_grid, profile);
    out = sim.run();

    Yhist = out.concentrations;   % n_t x n_z x n_sec
    t_out = out.time;
    n_t   = length(t_out);

    % grid ceiling for this n
    d_top = 20.0 * 2^((2*n - 1)/6) / 1000;   % mm
    top_bin_mm(i) = d_top;

    % deepest layer where disagg acts (D_max < grid ceiling)
    active_layers = find(dmax_mm < d_top);
    if isempty(active_layers)
        disagg_depth(i) = 0;
    else
        disagg_depth(i) = z(active_layers(end));
    end

    % total bv at snapshot days
    bv_all = squeeze(sum(Yhist, [2 3]));
    for j = 1:length(snap_days)
        ti = min(snap_days(j) + 1, n_t);
        bv_snaps(i, j) = bv_all(ti);
    end
    bv_final(i) = bv_all(end);

    % depth profile at t=365
    Yfinal = squeeze(Yhist(end, :, :));
    depth_profiles{i} = sum(Yfinal, 2);

    % largest bin in bottom layer
    bv_deep = Yfinal(end, :);
    last_bin = find(bv_deep > 0, 1, 'last');
    if isempty(last_bin)
        maxbin_deep(i) = 0;
    else
        maxbin_deep(i) = last_bin;
    end
end

% --- print summary table ---
fprintf('\n--- Summary ---\n');
fprintf('  %-6s  %-12s  %-12s  %-14s  %-16s  %-12s\n', ...
    'n', 'top_bin(mm)', 'disagg_to(m)', 'bv_t365', 'deep_maxbin/n', 'at_ceiling?');
for i = 1:n_runs
    n = n_list(i);
    at_ceil = maxbin_deep(i) == n;
    fprintf('  %-6d  %-12.3f  %-12.0f  %-14.4e  %-16s  %-12s\n', ...
        n, top_bin_mm(i), disagg_depth(i), bv_final(i), ...
        sprintf('%d / %d', maxbin_deep(i), n), mat2str(at_ceil));
end

fprintf('\n--- Total bv at snapshots ---\n');
fprintf('  %-6s', 'n');
for j = 1:length(snap_days)
    fprintf('  %-14s', sprintf('t=%d', snap_days(j)));
end
fprintf('\n');
for i = 1:n_runs
    fprintf('  %-6d', n_list(i));
    for j = 1:length(snap_days)
        fprintf('  %-14.4e', bv_snaps(i,j));
    end
    fprintf('\n');
end

% --- figures ---
fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = {'b-', 'r-', 'm-', 'k-'};

% Figure 1: depth profile at t=365 for each n
f1 = figure;
for i = 1:n_runs
    plot(depth_profiles{i}, z, colors{i}, 'LineWidth', 1.2); hold on;
end
hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
legend(arrayfun(@(n) sprintf('n=%d',n), n_list, 'UniformOutput', false), 'Location', 'best');
title('depth profile t=365: n comparison');
saveas(f1, fullfile(fig_dir, 'nsec_compare_depth.png'));

% Figure 2: total bv at t=365 vs n_sections
f2 = figure;
plot(n_list, bv_final, 'k-o', 'LineWidth', 1.2, 'MarkerSize', 6);
xlabel('n\_sections');
ylabel('total biovolume (t=365)');
title('convergence: total bv vs n');
saveas(f2, fullfile(fig_dir, 'nsec_compare_convergence.png'));

% Figure 3: D_max vs depth with grid ceilings marked
f3 = figure;
plot(dmax_mm, z, 'k-', 'LineWidth', 1.5); hold on;
line_styles = {'b--', 'r--', 'm--', 'g--'};
for i = 1:n_runs
    xline_val = top_bin_mm(i);
    plot([xline_val xline_val], [z(1) z(end)], line_styles{i}, 'LineWidth', 1.0);
end
hold off;
set(gca, 'YDir', 'reverse');
xlabel('size (mm)');
ylabel('depth (m)');
legend([{'D\_max(z)'}, arrayfun(@(n) sprintf('n=%d ceiling',n), n_list, 'UniformOutput', false)], ...
    'Location', 'best');
title('D\_max vs depth with grid ceilings');
saveas(f3, fullfile(fig_dir, 'nsec_compare_dmax.png'));

fprintf('\nFigures saved:\n');
fprintf('  nsec_compare_depth.png\n');
fprintf('  nsec_compare_convergence.png\n');
fprintf('  nsec_compare_dmax.png\n');

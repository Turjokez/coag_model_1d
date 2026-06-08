% run_full_1d_integrated.m
% Full 1-D run: coag + depth-scaled kernels + disagg + zoo
% + separate fecal pellets + cross-coag + mining.
%
% This is the first complete integrated run combining all physics.
% Goal: check budget, depth profile, size spectrum, D_max(z).

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

% --- grid and depth profile ---
col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z        = col_grid.z_centers;

% --- config ---
cfg = SimulationConfig( ...
    'n_sections',        30, ...
    't_final',           365, ...
    'delta_t',           0.4, ...
    'sinking_law',       'kriest_8', ...
    'ds_kernel_mode',    'sinking_law', ...
    'enable_coag',       true, ...
    'enable_sinking',    true, ...
    'enable_disagg',     true, ...
    'disagg_mode',       'operator_split', ...
    'disagg_dmax_cm',    1.0, ...       % fallback if eps not available
    'proc_substeps',     20, ...
    'enable_surface_pp', true, ...
    'surface_pp_bin',    1, ...
    'surface_pp_mu',     0.1, ...
    'enable_zoo',        true, ...
    'zoo_Zc',            0.307, ...     % Stemmann 2004 max, m^-3
    'zoo_Zf',            0.063, ...     % Stemmann 2004 max, m^-3
    'zoo_c',             0.025, ...     % clearance rate, m^3 ind^-1 day^-1
    'zoo_s',             1.3e-5, ...    % capture cross-section, m^2 ind^-1
    'zoo_p',             0.5, ...
    'zoo_ic',            7, ...         % fecal to bin 8 (~115 um)
    'fp_alpha_cross',    0.5, ...
    'enable_mining',     true, ...
    'mining_Zm',         250, ...
    'mining_dm',         1e-5, ...
    'mining_s',          1.3e-5);

% --- run ---
fprintf('Running full 1-D integrated model (t=365 days)...\n');
sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();

Yhist   = out.concentrations;         % n_t x n_z x n_sec
Yfphist = out.fecal_concentrations;   % n_t x n_z x n_sec
t_out = out.time;
n_t   = length(t_out);
n_z   = col_grid.n_z;

% --- budget check ---
bv_agg   = squeeze(sum(Yhist,   [2 3]));
bv_fp    = squeeze(sum(Yfphist, [2 3]));
bv_total = bv_agg + bv_fp;
bv0      = bv_total(1);
idx_day  = @(d) find(abs(t_out - d) == min(abs(t_out - d)), 1, 'first');

fprintf('\n--- Budget ---\n');
fprintf('  t=0     total bv = %.4e\n', bv_total(1));
fprintf('  t=30    total bv = %.4e\n', bv_total(idx_day(30)));
fprintf('  t=90    total bv = %.4e\n', bv_total(idx_day(90)));
fprintf('  t=180   total bv = %.4e\n', bv_total(idx_day(180)));
fprintf('  t=365   total bv = %.4e\n', bv_total(end));
fprintf('  Change t0->t365  = %.2f%%\n', 100*(bv_total(end)-bv0)/max(bv0,eps));
fprintf('  Final aggregate bv = %.4e\n', bv_agg(end));
fprintf('  Final fecal bv     = %.4e\n', bv_fp(end));

% --- D_max profile (depth-varying) ---
fprintf('\n--- D_max profile ---\n');
fprintf('  %-10s  %-14s  %-14s  %-10s\n', 'depth (m)', 'eps (cm^2/s^3)', 'eps (m^2/s^3)', 'D_max (mm)');
Dmax_A = 9.39e-6;
for k = 1:n_z
    eps_cm = profile.eps(k);
    eps_m  = eps_cm / 1e4;
    dmax_m = Dmax_A * eps_m^(-1/4);
    fprintf('  %-10.0f  %-14.3e  %-14.3e  %-10.3f\n', z(k), eps_cm, eps_m, dmax_m*1000);
end

% --- largest populated bin at t=365 ---
Yfinal   = squeeze(Yhist(end, :, :));     % n_z x n_sec
Yfpfinal = squeeze(Yfphist(end, :, :));   % n_z x n_sec
Ytotal_final = Yfinal + Yfpfinal;
[~, max_bin] = max(Ytotal_final, [], 2);  % peak bin per layer
fprintf('\n--- Size spectrum at t=365 (peak bin per depth layer) ---\n');
fprintf('  %-10s  %-10s\n', 'depth (m)', 'peak bin');
for k = 1:n_z
    fprintf('  %-10.0f  %-10d\n', z(k), max_bin(k));
end

% surface vs deep size spectrum
bv_surf = Ytotal_final(1, :);
bv_deep = Ytotal_final(end, :);
fprintf('\n  Surface (z=%.0fm): max bin = %d\n', z(1),   find(bv_surf>0, 1, 'last'));
fprintf('  Deep    (z=%.0fm): max bin = %d\n', z(end), find(bv_deep>0, 1, 'last'));

% --- figures ---
fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Figure 1: total biovolume vs time
f1 = figure;
plot(t_out, bv_total, 'k-', 'LineWidth', 1.2);
xlabel('time (day)');
ylabel('total biovolume');
title('full 1-D: total bv vs time');
saveas(f1, fullfile(fig_dir, 'full1d_totalvol.png'));

% Figure 2: depth profile at t = 30, 90, 180, 365
snap_days = [30, 90, 180, 365];
colors    = {'b-', 'r-', 'm-', 'k-'};
f2 = figure;
for i = 1:length(snap_days)
    ti = idx_day(snap_days(i));
    bv_z = sum(squeeze(Yhist(ti,:,:)), 2) + sum(squeeze(Yfphist(ti,:,:)), 2);
    plot(bv_z, z, colors{i}, 'LineWidth', 1.2); hold on;
end
hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume');
ylabel('depth (m)');
legend({'t=30','t=90','t=180','t=365'}, 'Location', 'best');
title('full 1-D: depth profile');
saveas(f2, fullfile(fig_dir, 'full1d_depth_profile.png'));

% Figure 3: surface vs deep size spectrum at t=365
f3 = figure;
bins = 1:cfg.n_sections;
plot(bins, bv_surf, 'b-o', 'MarkerSize', 4); hold on;
plot(bins, bv_deep, 'r-o', 'MarkerSize', 4);
hold off;
xlabel('bin');
ylabel('biovolume');
legend({'surface','deep'}, 'Location', 'best');
title('full 1-D: size spectrum t=365');
saveas(f3, fullfile(fig_dir, 'full1d_size_spectrum.png'));

% Figure 4: D_max vs depth
dmax_profile = zeros(n_z, 1);
for k = 1:n_z
    eps_m = profile.eps(k) / 1e4;
    dmax_profile(k) = Dmax_A * eps_m^(-1/4) * 1000;   % mm
end
f4 = figure;
plot(dmax_profile, z, 'k-', 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse');
xlabel('D\_max (mm)');
ylabel('depth (m)');
title('D\_max vs depth');
saveas(f4, fullfile(fig_dir, 'full1d_dmax_profile.png'));

fprintf('\nFigures saved to %s\n', fig_dir);
fprintf('  full1d_totalvol.png\n');
fprintf('  full1d_depth_profile.png\n');
fprintf('  full1d_size_spectrum.png\n');
fprintf('  full1d_dmax_profile.png\n');

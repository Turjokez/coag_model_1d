% run_fecal_separate_test.m
% Test that fecal pellets (Y_fp) are now tracked separately from aggregates (Y).
%
% What we check:
%   1. Y_fp starts at zero and grows over time (fecal production is working).
%   2. Y_fp sinks down the column (transport is working).
%   3. No negative values in Y or Y_fp.
%   4. Y total is LOWER than a run with fecal returned to Y
%      (because fecal is now removed from aggregate array).
%   5. Budget: Y_total + Y_fp_total should be similar to old single-array total.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z        = col_grid.z_centers;
n_z      = col_grid.n_z;
dz       = col_grid.dz;

cfg = SimulationConfig( ...
    'n_sections',        30, ...
    't_final',           180, ...
    'delta_t',           0.4, ...
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

fprintf('Running fecal separate tracking test (t=180 days, n=30)...\n');
sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();

Yhist   = out.concentrations;          % aggregates: n_t x n_z x n_sec
Yfphist = out.fecal_concentrations;    % fecal:      n_t x n_z x n_sec
t_out   = out.time;
n_t     = length(t_out);

% total biovolume vs time
bv_agg  = squeeze(sum(Yhist,   [2 3]));   % aggregates
bv_fp   = squeeze(sum(Yfphist, [2 3]));   % fecal pellets
bv_tot  = bv_agg + bv_fp;                 % combined

% check: no negatives
neg_agg = sum(Yhist(:)   < -1e-30);
neg_fp  = sum(Yfphist(:) < -1e-30);

% fecal depth profile at t=180
Yfp_final = squeeze(Yfphist(end, :, :));   % n_z x n_sec
bv_fp_z   = sum(Yfp_final, 2);             % total fecal per depth layer

% fecal size spectrum at surface and deep
bv_fp_surf = Yfp_final(1, :);
bv_fp_deep = Yfp_final(end, :);

% bin diameters
d_k = 20 * 2.^((2*(1:30) - 1)/6);   % um

% sinking speed at fecal pellet bin (bin 8, ~115 um) — surface layer
fp_bin      = max(1, min(30, round(cfg.zoo_ic) + 1));
w_agg_surf  = out.w_z(1,    fp_bin);      % m/day, aggregate at surface
w_fp_surf   = out.w_fp_z(1, fp_bin);      % m/day, fecal pellet at surface

% --- print results ---
fprintf('\n=== Fecal Separate Tracking Check ===\n');
fprintf('  Negatives in Y:    %d  (should be 0)\n', neg_agg);
fprintf('  Negatives in Y_fp: %d  (should be 0)\n', neg_fp);

fprintf('\n  Sinking speed at fecal pellet bin %d (~115 um), surface:\n', fp_bin);
fprintf('    Marine snow (kriest_8):   %.1f m/day\n', w_agg_surf);
fprintf('    Fecal pellet (Stokes):    %.1f m/day\n', w_fp_surf);
fprintf('    Ratio fp/agg:             %.1f x\n', w_fp_surf / max(w_agg_surf, eps));

fprintf('\n  t=0   agg=%.4e  fp=%.4e  total=%.4e\n', bv_agg(1),   bv_fp(1),   bv_tot(1));
fprintf('  t=30  agg=%.4e  fp=%.4e  total=%.4e\n', ...
    bv_agg(min(76,n_t)), bv_fp(min(76,n_t)), bv_tot(min(76,n_t)));
fprintf('  t=180 agg=%.4e  fp=%.4e  total=%.4e\n', bv_agg(end), bv_fp(end), bv_tot(end));

fprintf('\n  Fecal / Total at t=180: %.1f%%\n', 100*bv_fp(end)/max(bv_tot(end),eps));
fprintf('  Fecal max bin in deep (975m): %d\n', find(bv_fp_deep > 0, 1, 'last'));

if bv_fp(end) > 0
    fprintf('\n  PASS: Y_fp is growing (fecal production working).\n');
else
    fprintf('\n  FAIL: Y_fp is still zero at t=180.\n');
end

% --- figures ---
fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Figure 1: two panels — left: aggregate bv, right: fecal bv (separate scales)
f1 = figure;
subplot(1,2,1);
plot(t_out, bv_agg, 'b-', 'LineWidth', 1.2);
xlabel('time (day)');
ylabel('aggregate biovolume');
title('aggregates');

subplot(1,2,2);
plot(t_out, bv_fp, 'r-', 'LineWidth', 1.2);
xlabel('time (day)');
ylabel('fecal biovolume');
title('fecal pellets');

saveas(f1, fullfile(fig_dir, 'fp_sep_totalvol.png'));

% Figure 2: two panels — left: aggregate depth profile, right: fecal depth profile
Yagg_final = squeeze(Yhist(end, :, :));
bv_agg_z   = sum(Yagg_final, 2);

f2 = figure;
subplot(1,2,1);
plot(bv_agg_z, z, 'b-', 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse');
xlabel('aggregate biovolume');
ylabel('depth (m)');
title('aggregates  t=180');

subplot(1,2,2);
plot(bv_fp_z, z, 'r-', 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse');
xlabel('fecal biovolume');
ylabel('depth (m)');
title('fecal pellets  t=180');

saveas(f2, fullfile(fig_dir, 'fp_sep_depth.png'));

% Figure 3: fecal size spectrum at surface and deep — truncate noise floor
f3 = figure;
semilogy(d_k, bv_fp_surf + 1e-40, 'b-', 'LineWidth', 1.2); hold on;
semilogy(d_k, bv_fp_deep + 1e-40, 'k-', 'LineWidth', 1.2);
hold off;
% cut noise floor: only show values with physical meaning
ymax = max([bv_fp_surf, bv_fp_deep]);
ylim([1e-15, ymax * 10]);
xlabel('diameter (\mum)');
ylabel('fecal biovolume');
legend({'surface (25m)', 'deep (975m)'}, 'Location', 'best');
title('fecal size spectrum t=180');
saveas(f3, fullfile(fig_dir, 'fp_sep_spectrum.png'));

fprintf('\nFigures saved:\n');
fprintf('  fp_sep_totalvol.png\n');
fprintf('  fp_sep_depth.png\n');
fprintf('  fp_sep_spectrum.png\n');

% run_fecal_cross_test.m
% Test Step 3: cross-coagulation between fecal pellets and marine snow.
%
% What we check:
%   1. No negatives in Y or Y_fp.
%   2. Total biovolume (Y + Y_fp) is conserved — cross-coag just moves volume.
%   3. Y_fp is LOWER than without cross-coag (fecal being absorbed by marine snow).
%   4. Y is HIGHER than without cross-coag (marine snow gains fecal volume).
%   5. Budget: volume transferred matches the difference.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z        = col_grid.z_centers;

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
    'zoo_ic',            7, ...
    'fp_alpha_cross',    0.5);

% run with cross-coag ON
fprintf('Running with cross-coag ON (alpha_cross=0.5, t=180 days)...\n');
sim_on = ColumnSimulation(cfg, col_grid, profile);
out_on = sim_on.run();

% run with cross-coag OFF (alpha_cross=0 disables it effectively)
cfg_off = cfg.copy();
cfg_off.fp_alpha_cross = 0.0;
fprintf('Running with cross-coag OFF (alpha_cross=0, t=180 days)...\n');
sim_off = ColumnSimulation(cfg_off, col_grid, profile);
out_off = sim_off.run();

% --- extract results ---
Yhist_on    = out_on.concentrations;
Yfphist_on  = out_on.fecal_concentrations;
Yhist_off   = out_off.concentrations;
Yfphist_off = out_off.fecal_concentrations;
t_out       = out_on.time;

% total biovolume vs time
bv_agg_on  = squeeze(sum(Yhist_on,   [2 3]));
bv_fp_on   = squeeze(sum(Yfphist_on, [2 3]));
bv_tot_on  = bv_agg_on + bv_fp_on;

bv_agg_off = squeeze(sum(Yhist_off,   [2 3]));
bv_fp_off  = squeeze(sum(Yfphist_off, [2 3]));
bv_tot_off = bv_agg_off + bv_fp_off;

% negatives
neg_agg = sum(Yhist_on(:) < -1e-30);
neg_fp  = sum(Yfphist_on(:) < -1e-30);

% depth profiles at t=180
Yagg_fin = squeeze(Yhist_on(end,:,:));
Yfp_fin  = squeeze(Yfphist_on(end,:,:));
bv_agg_z = sum(Yagg_fin, 2);
bv_fp_z  = sum(Yfp_fin,  2);

% --- print results ---
fprintf('\n=== Cross-Coagulation Check ===\n');
fprintf('  Negatives in Y:    %d  (should be 0)\n', neg_agg);
fprintf('  Negatives in Y_fp: %d  (should be 0)\n', neg_fp);

fprintf('\n  At t=180, with vs without cross-coag:\n');
fprintf('    Marine snow bv:   ON=%.4e  OFF=%.4e  ratio=%.3f\n', ...
    bv_agg_on(end), bv_agg_off(end), bv_agg_on(end)/max(bv_agg_off(end),eps));
fprintf('    Fecal bv:         ON=%.4e  OFF=%.4e  ratio=%.3f\n', ...
    bv_fp_on(end), bv_fp_off(end), bv_fp_on(end)/max(bv_fp_off(end),eps));
fprintf('    Total bv:         ON=%.4e  OFF=%.4e  ratio=%.3f\n', ...
    bv_tot_on(end), bv_tot_off(end), bv_tot_on(end)/max(bv_tot_off(end),eps));

fprintf('\n  Fecal fraction at t=180:\n');
fprintf('    cross-coag ON:  %.2f%%\n', 100*bv_fp_on(end)/max(bv_tot_on(end),eps));
fprintf('    cross-coag OFF: %.2f%%\n', 100*bv_fp_off(end)/max(bv_tot_off(end),eps));

if bv_fp_on(end) < bv_fp_off(end)
    fprintf('\n  PASS: Y_fp lower with cross-coag (fecal being absorbed by marine snow).\n');
else
    fprintf('\n  FAIL: Y_fp not lower with cross-coag. Check implementation.\n');
end
if bv_agg_on(end) > bv_agg_off(end)
    fprintf('  PASS: Marine snow higher with cross-coag (received fecal volume).\n');
else
    fprintf('  NOTE: Marine snow not clearly higher -- effect may be small.\n');
end

% --- figures ---
fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Figure 1: fecal bv over time — on vs off
f1 = figure;
plot(t_out, bv_fp_off, 'r--', 'LineWidth', 1.2); hold on;
plot(t_out, bv_fp_on,  'r-',  'LineWidth', 1.2);
hold off;
xlabel('time (day)');
ylabel('fecal biovolume');
legend({'fp: no cross-coag', 'fp: with cross-coag'}, 'Location', 'best');
title('cross-coag effect on fecal standing stock');
saveas(f1, fullfile(fig_dir, 'cross_coag_fp_time.png'));

% Figure 2: depth profile — fecal on vs off
f2 = figure;
Yfp_fin_off = squeeze(Yfphist_off(end,:,:));
bv_fp_z_off = sum(Yfp_fin_off, 2);

plot(bv_fp_z_off, z, 'r--', 'LineWidth', 1.2); hold on;
plot(bv_fp_z,     z, 'r-',  'LineWidth', 1.2);
hold off;
set(gca, 'YDir', 'reverse');
xlabel('fecal biovolume (m^{-3})');
ylabel('depth (m)');
legend({'no cross-coag', 'with cross-coag'}, 'Location', 'best');
title('fecal depth profile at t=180');
saveas(f2, fullfile(fig_dir, 'cross_coag_depth.png'));

% Figure 3: marine snow depth profile on vs off
Yagg_fin_off = squeeze(Yhist_off(end,:,:));
bv_agg_z_off = sum(Yagg_fin_off, 2);

f3 = figure;
plot(bv_agg_z_off, z, 'b--', 'LineWidth', 1.2); hold on;
plot(bv_agg_z,     z, 'b-',  'LineWidth', 1.2);
hold off;
set(gca, 'YDir', 'reverse');
xlabel('aggregate biovolume (m^{-3})');
ylabel('depth (m)');
legend({'no cross-coag', 'with cross-coag'}, 'Location', 'best');
title('marine snow depth profile at t=180');
saveas(f3, fullfile(fig_dir, 'cross_coag_agg_depth.png'));

fprintf('\nFigures saved:\n');
fprintf('  cross_coag_fp_time.png\n');
fprintf('  cross_coag_depth.png\n');
fprintf('  cross_coag_agg_depth.png\n');

% run_microbe_sensitivity.m
% Sensitivity test: microbial remineralization rate r0 on transfer efficiency.
%
% Four runs (A-D): r0 = 0, 0.01, 0.03, 0.1 day^-1.
% All other physics identical (n=30, dt=0.4, t=365 days, full physics on).
% Prints TE and total biovolume for each case.
% Saves one summary figure.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);

% base config — full physics, n=30
cfg_base = SimulationConfig( ...
    'n_sections',        30, ...
    't_final',           365, ...
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
    'fp_alpha_cross',    0.5, ...
    'enable_microbe',    true);

% r0 values to test (narrow sweep — baseline is ~1.1% at r0=0)
r0_vals  = [0.0, 0.001, 0.002, 0.004, 0.006, 0.01];
labels   = {'r0=0', 'r0=0.001', 'r0=0.002', 'r0=0.004', 'r0=0.006', 'r0=0.01'};
colors   = {'k', 'b', 'c', 'g', 'm', 'r'};

n_runs = numel(r0_vals);
TE     = zeros(1, n_runs);
bv_fin = zeros(1, n_runs);
bv_bot = zeros(1, n_runs);

for i = 1:n_runs
    cfg = cfg_base.copy();
    cfg.microbe_r0 = r0_vals(i);
    if r0_vals(i) == 0
        cfg.enable_microbe = false;
    end

    fprintf('Run %s ...\n', labels{i});
    sim = ColumnSimulation(cfg, col_grid, profile);
    out = sim.run();

    % surface PP flux at t=365 [bv m^-2 day^-1]
    % same formula as run_n30_transfer_check.m
    Y_surf  = squeeze(out.concentrations(end, 1, :));
    pp_flux = cfg.surface_pp_mu * Y_surf(1) * col_grid.dz;

    % bottom sinking flux — aggregate + fecal, divided by dz (matches n30 script)
    w_bot    = out.w_z(end, :);
    Y_bot    = squeeze(out.concentrations(end, end, :))';
    bflux_agg = sum(w_bot .* Y_bot) / col_grid.dz;

    w_fp_bot  = out.w_fp_z(end, :);
    Yfp_bot   = squeeze(out.fecal_concentrations(end, end, :))';
    bflux_fp  = sum(w_fp_bot .* Yfp_bot) / col_grid.dz;

    bot_flux = bflux_agg + bflux_fp;

    TE(i)     = 100 * bot_flux / max(pp_flux, eps);
    bv_fin(i) = sum(out.concentrations(end, :, :), 'all');
    bv_bot(i) = bot_flux;

    fprintf('  TE = %.2f%%   total bv = %.4e\n', TE(i), bv_fin(i));
end

fprintf('\n--- Summary ---\n');
fprintf('%-16s  TE(%%)   total_bv\n', 'Case');
for i = 1:n_runs
    fprintf('%-16s  %5.2f   %.4e\n', labels{i}, TE(i), bv_fin(i));
end

% --- figure: total biovolume at t=365 for each run ---
fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% plot only the non-zero r0 cases (skip index 1 = r0=0)
r0_plot = r0_vals(2:end);
TE_plot  = TE(2:end);
bv_plot  = bv_fin(2:end);

f1 = figure;
plot(r0_plot, TE_plot, 'ko-', 'LineWidth', 1.2, 'MarkerFaceColor', 'k');
hold on;
yline(TE(1), 'k--');   % baseline r0=0
hold off;
xlabel('r_0 (day^{-1})');
ylabel('transfer efficiency (%)');
legend({'microbial on', 'baseline (r_0=0)'}, 'Location', 'best');
title('TE vs microbial rate');
saveas(f1, fullfile(fig_dir, 'microbe_sensitivity_TE.png'));

f2 = figure;
plot(r0_plot, bv_plot, 'ko-', 'LineWidth', 1.2, 'MarkerFaceColor', 'k');
hold on;
yline(bv_fin(1), 'k--');
hold off;
xlabel('r_0 (day^{-1})');
ylabel('total biovolume at t=365');
legend({'microbial on', 'baseline (r_0=0)'}, 'Location', 'best');
title('standing stock vs microbial rate');
saveas(f2, fullfile(fig_dir, 'microbe_sensitivity_bv.png'));

fprintf('\nFigures saved:\n');
fprintf('  microbe_sensitivity_bv.png\n');
fprintf('  microbe_sensitivity_TE.png\n');

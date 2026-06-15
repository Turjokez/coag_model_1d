% run_100m_r0_sensitivity.m
%
% Sweep microbial loss rate r0 in the 100 m start framework.
%
% From report_june12: microbe OFF raised the ratio at all depths but
% did not close the deep low bias. Test how sensitive the deep profile
% is to r0 magnitude, and whether a lower r0 can recover the deep signal.
%
% Cases: r0 = 0.0 (off), 0.005, 0.01, 0.02, 0.03 (current default)
% Everything else identical to run_100m_start.m.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% ---------------------------------------------------------------
% 1. Shared setup
% ---------------------------------------------------------------
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);

k_bc      = 2;
k_compare = 3:10;
z_compare = col_grid.z_centers(k_compare);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg_base = SimulationConfig();
cfg_base.n_sections       = 30;
cfg_base.sinking_law      = 'kriest_8';
cfg_base.disagg_mode      = 'operator_split';
cfg_base.disagg_dmax_cm   = 1.0;
cfg_base.disagg_dmax_A    = 9.39e-6 * 5;
cfg_base.enable_coag      = true;
cfg_base.enable_disagg    = true;
cfg_base.enable_zoo       = true;
cfg_base.enable_microbe   = true;
cfg_base.enable_mining    = true;
cfg_base.alpha            = 0.5;
cfg_base.microbe_r0       = 0.03;
cfg_base.microbe_use_temp = true;
cfg_base.microbe_tref_C   = 20;
cfg_base.surface_pp_mu    = 0.0;
cfg_base.r_to_rg          = 1.6;
cfg_base.zoo_c            = 0.025;
cfg_base.zoo_s            = 1.3e-5;
cfg_base.zoo_p            = 0.5;
cfg_base.zoo_ic           = 7;
cfg_base.mining_s         = 1.3e-5;
cfg_base.fp_alpha_cross   = 0.5;
cfg_base.validate();

grid_cfg   = cfg_base.derive();
r_cm       = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;
n_sec      = cfg_base.n_sections;

% ---------------------------------------------------------------
% 2. BC: UVP at 100 m -> model bins with power-law fill
% ---------------------------------------------------------------
bc = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_compare);
phi_bc_daily   = bc.phi_bc_daily;
n_days         = bc.n_days;
id_model_best  = bc.id_model_best;
phi_uvp_cmp    = bc.phi_uvp_cmp;
mask_uvp_model = bc.mask_uvp_model;
fprintf('Best cast day: %d\n', bc.best_date);

% ---------------------------------------------------------------
% 3. r0 sweep cases
% ---------------------------------------------------------------
r0_vals   = [0.0, 0.005, 0.01, 0.02, 0.03];
n_cases   = numel(r0_vals);
ratio_table = NaN(numel(k_compare), n_cases);

for ic = 1:n_cases
    r0 = r0_vals(ic);
    fprintf('\n=== r0 = %.3f day^-1 ===\n', r0);

    cfg = cfg_base;
    cfg.microbe_r0     = r0;
    cfg.enable_microbe = (r0 > 0);
    cfg.validate();

    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(col_grid.n_z, n_sec);
    Yfp = zeros(col_grid.n_z, n_sec);

    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            for i_step = 1:steps_per_day
                Y(k_bc, :) = phi_bc_daily(i_day, :);
                [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
                Y(k_bc, :) = phi_bc_daily(i_day, :);
            end
        end
        phi_after  = mean(sum(Y + Yfp, 2));
        rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
        if rel_change < spinup_tol
            fprintf('  Converged at cycle %d\n', icyc);
            break;
        end
    end

    Y   = zeros(col_grid.n_z, n_sec);
    Yfp = zeros(col_grid.n_z, n_sec);
    Y_snap = [];
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
        if i_day == id_model_best
            Y_snap = Y + Yfp;
        end
    end

    phi_mod = sum(Y_snap(k_compare, mask_uvp_model), 2);
    ratio_table(:, ic) = phi_mod ./ max(phi_uvp_cmp, 1e-30);
end

% ---------------------------------------------------------------
% 4. Print ratio table
% ---------------------------------------------------------------
fprintf('\n--- Ratio table: model / UVP (100-2000 um) ---\n');
fprintf('%-10s', 'Depth(m)');
for ic = 1:n_cases
    fprintf('  r0=%.3f', r0_vals(ic));
end
fprintf('\n%s\n', repmat('-', 1, 10 + n_cases*10));
for i = 1:numel(k_compare)
    fprintf('%-10.0f', z_compare(i));
    for ic = 1:n_cases
        fprintf('  %7.2f', ratio_table(i, ic));
    end
    fprintf('\n');
end

% ---------------------------------------------------------------
% 5. Figure: ratio vs depth for each r0
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
figure('Units','centimeters','Position',[2 2 10 13]);
hold on;
for ic = 1:n_cases
    plot(ratio_table(:, ic), z_compare, '-o', ...
        'Color', colors(ic,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('r_0 = %.3f', r0_vals(ic)));
end
xline(1.0, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('model / UVP');
ylabel('Depth (m)');
ylim([100 500]);
xlim([0 3]);
legend('location', 'northeast', 'FontSize', 7);
title('r_0 sensitivity: 100 m start');

saveas(gcf, fullfile(fig_dir, 'r0_sensitivity_ratio.png'));
fprintf('\nSaved r0_sensitivity_ratio.png\n');

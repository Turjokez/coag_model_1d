% run_100m_alpha_sensitivity.m
%
% Sweep coagulation efficiency alpha in the 100 m start framework.
%
% From report_june12 Section 10: adding coagulation causes the biggest
% single drop in the ratio. At alpha=0.5 the model loses too much mass
% from the UVP window. Test whether lower alpha recovers the deep signal.
%
% Cases: alpha = 0.05, 0.10, 0.20, 0.30, 0.50 (current default)
% All other settings identical to run_100m_start (full physics, r0=0.03).

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
% 2. BC at 100 m with power-law fill below 100 um
% ---------------------------------------------------------------
bc = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_compare);
phi_bc_daily   = bc.phi_bc_daily;
n_days         = bc.n_days;
id_model_best  = bc.id_model_best;
phi_uvp_cmp    = bc.phi_uvp_cmp;
mask_uvp_model = bc.mask_uvp_model;
fprintf('Best cast day: %d\n', bc.best_date);

% ---------------------------------------------------------------
% 3. Alpha sweep
% ---------------------------------------------------------------
alpha_vals  = [0.05, 0.10, 0.20, 0.30, 0.50];
n_cases     = numel(alpha_vals);
ratio_table = NaN(numel(k_compare), n_cases);
Y_snaps     = cell(n_cases, 1);

for ic = 1:n_cases
    al = alpha_vals(ic);
    fprintf('\n=== alpha = %.2f ===\n', al);

    cfg       = cfg_base;
    cfg.alpha = al;
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
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
        if i_day == id_model_best
            Y_snaps{ic} = Y + Yfp;
        end
    end

    phi_mod = sum(Y_snaps{ic}(k_compare, mask_uvp_model), 2);
    ratio_table(:, ic) = phi_mod ./ max(phi_uvp_cmp, 1e-30);
end

% ---------------------------------------------------------------
% 4. Print ratio table
% ---------------------------------------------------------------
fprintf('\n--- Alpha sensitivity: model (UVP range) / UVP ---\n');
fprintf('%-10s', 'Depth(m)');
for ic = 1:n_cases
    fprintf('  alpha=%.2f', alpha_vals(ic));
end
fprintf('\n%s\n', repmat('-', 1, 10 + n_cases*12));
for i = 1:numel(k_compare)
    fprintf('%-10.0f', z_compare(i));
    for ic = 1:n_cases
        fprintf('  %9.2f', ratio_table(i, ic));
    end
    fprintf('\n');
end

% ---------------------------------------------------------------
% 5. Figure: ratio vs depth for each alpha
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
figure('Units','centimeters','Position',[2 2 11 13]);
hold on;
for ic = 1:n_cases
    plot(ratio_table(:, ic), z_compare, '-o', ...
        'Color', colors(ic,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('\\alpha = %.2f', alpha_vals(ic)));
end
xline(1.0, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('model / UVP');
ylabel('Depth (m)');
ylim([100 500]);
xlim([0 2.5]);
legend('location', 'northeast', 'FontSize', 7);
title('\alpha sensitivity: 100 m start');

saveas(gcf, fullfile(fig_dir, 'alpha_sensitivity_ratio.png'));
fprintf('\nSaved alpha_sensitivity_ratio.png\n');

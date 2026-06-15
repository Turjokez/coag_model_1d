% run_100m_da_sweep.m
%
% Test whether the disagg Da multiplier controls the remaining deep bias.
%
% From report_june12 Section 12: best (alpha, r0) = (0.10, 0) fixes
% 125-325 m but leaves ratio 0.45-0.82 at 375-475 m. The 2D grid
% score surface is flat at the low-alpha corner, meaning alpha and r0
% alone cannot close the deep gap.
%
% Hypothesis: disagg_dmax_A = 9.39e-6 * 5 (Da x5) is too aggressive
% at depth. It fragments large particles into sub-100 um bins faster
% than coagulation can rebuild them. Lowering Da should keep more mass
% in the UVP-visible range at depth.
%
% Sweep: Da multiplier = [1, 2, 3, 5, 8]
% Fixed: alpha = 0.10, r0 = 0.0 (best pair from Section 12)

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

Da_base = 9.39e-6;   % Parker et al. base value [m]

cfg_base = SimulationConfig();
cfg_base.n_sections       = 30;
cfg_base.sinking_law      = 'kriest_8';
cfg_base.disagg_mode      = 'operator_split';
cfg_base.disagg_dmax_cm   = 1.0;
cfg_base.disagg_dmax_A    = Da_base * 5;   % current default
cfg_base.enable_coag      = true;
cfg_base.enable_disagg    = true;
cfg_base.enable_zoo       = true;
cfg_base.enable_microbe   = false;   % r0 = 0 (best from 2D grid)
cfg_base.enable_mining    = true;
cfg_base.alpha            = 0.10;   % best from 2D grid
cfg_base.microbe_r0       = 0.0;
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
% 2. BC at 100 m with power-law fill
% ---------------------------------------------------------------
bc = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_compare);
phi_bc_daily   = bc.phi_bc_daily;
n_days         = bc.n_days;
id_model_best  = bc.id_model_best;
phi_uvp_cmp    = bc.phi_uvp_cmp;
mask_uvp_model = bc.mask_uvp_model;
fprintf('Best cast day: %d\n', bc.best_date);

% ---------------------------------------------------------------
% 3. Da multiplier sweep
% ---------------------------------------------------------------
da_mults    = [1, 2, 3, 5, 8];
n_cases     = numel(da_mults);
ratio_table = NaN(numel(k_compare), n_cases);
scores      = NaN(1, n_cases);
Y_snaps     = cell(n_cases, 1);

for ic = 1:n_cases
    mult = da_mults(ic);
    fprintf('\n=== Da x%d ===\n', mult);

    cfg = cfg_base;
    cfg.disagg_dmax_A = Da_base * mult;
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
    scores(ic) = mean((ratio_table(:, ic) - 1).^2);
end

% ---------------------------------------------------------------
% 4. Print ratio table and scores
% ---------------------------------------------------------------
fprintf('\n--- Da sweep: model / UVP  (alpha=0.10, r0=0) ---\n');
fprintf('%-10s', 'Depth(m)');
for ic = 1:n_cases
    fprintf('  Da x%-4d', da_mults(ic));
end
fprintf('\n%s\n', repmat('-', 1, 10 + n_cases*10));
for i = 1:numel(k_compare)
    fprintf('%-10.0f', z_compare(i));
    for ic = 1:n_cases
        fprintf('  %7.2f', ratio_table(i, ic));
    end
    fprintf('\n');
end
fprintf('\nScores:');
for ic = 1:n_cases
    fprintf('  Da x%d: %.4f', da_mults(ic), scores(ic));
end
[~, ibest_da] = min(scores);
fprintf('\nBest Da multiplier: x%d\n', da_mults(ibest_da));

% ---------------------------------------------------------------
% 5. Figure: ratio vs depth
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
figure('Units','centimeters','Position',[2 2 11 13]);
hold on;
for ic = 1:n_cases
    plot(ratio_table(:, ic), z_compare, '-o', ...
        'Color', colors(ic,:), 'LineWidth', 1.5, ...
        'DisplayName', sprintf('D_a \\times%d', da_mults(ic)));
end
xline(1.0, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('model / UVP');
ylabel('Depth (m)');
ylim([100 500]);
xlim([0 2.5]);
legend('location', 'northeast', 'FontSize', 7);
title('D_a sweep: \alpha=0.10, r_0=0');

saveas(gcf, fullfile(fig_dir, 'da_sweep_ratio.png'));
fprintf('\nSaved da_sweep_ratio.png\n');

% run_alldredge_disagg_sweep.m
%
% Sweep disagg coefficient C0 across Parker-to-Alldredge range.
%
% Background:
%   Current model uses Parker et al. (1972): C0 = 9.39e-6 m
%   (derived from wastewater activated sludge).
%   Adrian suggests Alldredge marine-snow values which are larger
%   (weaker disagg per unit turbulence) because marine snow is
%   less dense and more fragile than wastewater flocs.
%
%   D_max = C0 * eps^(-1/4)
%
%   Values tested (all in metres):
%     9.39e-6   Parker original
%     1.88e-5   Parker x2
%     4.70e-5   Parker x5  <- current best from Da sweep
%     9.39e-5   Parker x10
%     2.35e-4   Parker x25
%     4.70e-4   Parker x50 (upper Alldredge-scale estimate)
%
% Fixed: alpha = 0.10, r0 = 0 (best from 2D grid, June 12)

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

Da_parker = 9.39e-6;   % Parker base [m]

cfg_base = SimulationConfig();
cfg_base.n_sections       = 30;
cfg_base.sinking_law      = 'kriest_8';
cfg_base.disagg_mode      = 'operator_split';
cfg_base.disagg_dmax_cm   = 1.0;
cfg_base.disagg_dmax_A    = Da_parker * 5;   % starting point
cfg_base.enable_coag      = true;
cfg_base.enable_disagg    = true;
cfg_base.enable_zoo       = true;
cfg_base.enable_microbe   = false;
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

n_sec = cfg_base.n_sections;

% ---------------------------------------------------------------
% 2. BC at 100 m
% ---------------------------------------------------------------
bc = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_compare);
phi_bc_daily   = bc.phi_bc_daily;
n_days         = bc.n_days;
id_model_best  = bc.id_model_best;
phi_uvp_cmp    = bc.phi_uvp_cmp;
mask_uvp_model = bc.mask_uvp_model;
fprintf('Best cast day: %d\n', bc.best_date);

% ---------------------------------------------------------------
% 3. C0 sweep
% ---------------------------------------------------------------
% Parker x multipliers: 1, 2, 5 (current), 10, 25, 50
c0_vals  = Da_parker * [1, 2, 5, 10, 25, 50];
c0_names = {'Parker x1', 'Parker x2', 'Parker x5 (current)', ...
             'Parker x10', 'Parker x25', 'Parker x50'};
n_cases  = numel(c0_vals);

ratio_table = NaN(numel(k_compare), n_cases);
scores      = NaN(1, n_cases);
Y_snaps     = cell(n_cases, 1);

for ic = 1:n_cases
    fprintf('\n=== %s  (C0 = %.2e m) ===\n', c0_names{ic}, c0_vals(ic));

    cfg = cfg_base;
    cfg.disagg_dmax_A = c0_vals(ic);
    cfg.validate();

    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(col_grid.n_z, n_sec);
    Yfp = zeros(col_grid.n_z, n_sec);

    % spinup
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

    % final run
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
% 4. Print results
% ---------------------------------------------------------------
fprintf('\n--- Alldredge sweep: model / UVP  (alpha=0.10, r0=0) ---\n');
fprintf('%-10s', 'Depth(m)');
for ic = 1:n_cases
    fprintf('  %8s', sprintf('x%d', round(c0_vals(ic)/Da_parker)));
end
fprintf('\n%s\n', repmat('-', 1, 10 + n_cases*10));
for i = 1:numel(k_compare)
    fprintf('%-10.0f', z_compare(i));
    for ic = 1:n_cases
        fprintf('  %8.2f', ratio_table(i, ic));
    end
    fprintf('\n');
end
fprintf('\nScore (lower=better):\n');
for ic = 1:n_cases
    fprintf('  %-25s  %.4f\n', c0_names{ic}, scores(ic));
end
[~, ibest] = min(scores);
fprintf('\nBest: %s\n', c0_names{ibest});

% ---------------------------------------------------------------
% 5. Figure
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
figure('Units','centimeters','Position',[2 2 11 13]);
hold on;
for ic = 1:n_cases
    plot(ratio_table(:, ic), z_compare, '-o', ...
        'Color', colors(ic,:), 'LineWidth', 1.5, ...
        'DisplayName', c0_names{ic});
end
xline(1.0, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('model / UVP');
ylabel('Depth (m)');
ylim([100 500]);
xlim([0 2.5]);
legend('location', 'northeast', 'FontSize', 7);
title('Alldredge C_0 sweep: \alpha=0.10, r_0=0');

saveas(gcf, fullfile(fig_dir, 'alldredge_disagg_sweep.png'));
fprintf('\nSaved alldredge_disagg_sweep.png\n');

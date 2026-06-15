% run_100m_2d_grid.m
%
% Small 2D grid for alpha and r0 in the 100 m start setup.
%
% Goal:
%   find the best (alpha, r0) pair for the 125-475 m UVP profile.
%
% Score:
%   mean squared relative error = mean((model/UVP - 1).^2)
%
% Outputs:
%   - prints score table
%   - prints best pair
%   - saves 2D score figure
%   - saves best model vs UVP profile figure

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

alpha_vals = [0.05, 0.07, 0.10, 0.15];
r0_vals    = [0.0, 0.005, 0.01, 0.02];

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
mask_uvp_model = d_model_um >= 100 & d_model_um < 2000;

% ---------------------------------------------------------------
% 2. Build 100 m BC from UVP with power-law fill
% ---------------------------------------------------------------
bc = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, k_compare);
phi_bc_daily  = bc.phi_bc_daily;
n_days        = bc.n_days;
id_model_best = bc.id_model_best;
phi_uvp_cmp   = bc.phi_uvp_cmp;
fprintf('Best cast day: %d\n', bc.best_date);

% ---------------------------------------------------------------
% 3. 2D grid
% ---------------------------------------------------------------
score_mat   = NaN(numel(alpha_vals), numel(r0_vals));
ratio_store = cell(numel(alpha_vals), numel(r0_vals));
best_score  = inf;
best_i = 1; best_j = 1;
best_phi = [];

for ia = 1:numel(alpha_vals)
    for ir = 1:numel(r0_vals)
        al = alpha_vals(ia);
        r0 = r0_vals(ir);
        fprintf('\n=== alpha = %.2f, r0 = %.3f ===\n', al, r0);

        cfg = cfg_base;
        cfg.alpha            = al;
        cfg.microbe_r0       = r0;
        cfg.enable_microbe   = (r0 > 0);
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
        ratio    = phi_mod ./ max(phi_uvp_cmp, 1e-30);
        score    = mean((ratio - 1).^2);

        ratio_store{ia, ir} = ratio;
        score_mat(ia, ir)   = score;

        fprintf('  score = %.4f\n', score);

        if score < best_score
            best_score = score;
            best_i     = ia;
            best_j     = ir;
            best_phi   = phi_mod;
        end
    end
end

% ---------------------------------------------------------------
% 4. Print summary
% ---------------------------------------------------------------
fprintf('\n--- 2D score table: mean((model/UVP - 1)^2) ---\n');
fprintf('%-12s', 'alpha\\r0');
for ir = 1:numel(r0_vals)
    fprintf('  %8.3f', r0_vals(ir));
end
fprintf('\n%s\n', repmat('-', 1, 12 + numel(r0_vals)*10));
for ia = 1:numel(alpha_vals)
    fprintf('%-12.2f', alpha_vals(ia));
    for ir = 1:numel(r0_vals)
        fprintf('  %8.4f', score_mat(ia, ir));
    end
    fprintf('\n');
end

fprintf('\nBest pair: alpha = %.2f, r0 = %.3f, score = %.4f\n', ...
    alpha_vals(best_i), r0_vals(best_j), best_score);

best_ratio = ratio_store{best_i, best_j};
fprintf('\nBest-pair ratio profile:\n');
for i = 1:numel(k_compare)
    fprintf('  z=%3.0f m : %.2f\n', z_compare(i), best_ratio(i));
end

% ---------------------------------------------------------------
% 5. Figures
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Figure 1: 2D score grid
figure('Units','centimeters','Position',[2 2 11 9]);
imagesc(r0_vals, alpha_vals, score_mat);
set(gca, 'YDir', 'normal');
hold on;
plot(r0_vals(best_j), alpha_vals(best_i), 'wx', 'MarkerSize', 10, 'LineWidth', 2);
xlabel('r_0 [day^{-1}]');
ylabel('\alpha');
title('100 m start: 2D score grid');
colorbar;
saveas(gcf, fullfile(fig_dir, 'grid2d_score.png'));

% Figure 2: best profile
figure('Units','centimeters','Position',[2 2 10 13]);
plot(best_phi, z_compare, 'b-o', 'LineWidth', 1.5, ...
    'DisplayName', sprintf('model: \\alpha=%.2f, r_0=%.3f', ...
    alpha_vals(best_i), r0_vals(best_j)));
hold on;
plot(phi_uvp_cmp, z_compare, 'k--', 'LineWidth', 1.5, 'DisplayName', 'UVP');
set(gca, 'YDir', 'reverse');
xlabel('\phi [ppmV]');
ylabel('Depth (m)');
ylim([100 500]);
legend('Location', 'southeast', 'FontSize', 7);
title('100 m start: best model vs UVP');
saveas(gcf, fullfile(fig_dir, 'grid2d_best_profile.png'));

fprintf('\nSaved grid2d_score.png\n');
fprintf('Saved grid2d_best_profile.png\n');

% Short note:
% 1. Build 100 m BC from UVP.
% 2. Run small alpha x r0 grid with same spinup.
% 3. Score each case by mean squared relative error.
% 4. Save the best pair and one clean profile figure.

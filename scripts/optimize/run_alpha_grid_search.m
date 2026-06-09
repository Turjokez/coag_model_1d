% run_alpha_grid_search.m
% Grid search over stickiness (alpha) against real EXPORTS UVP data.
%
% Loss = sum of squared log10 differences between model phi(z) and UVP phi(z).
% Log-space weights all depth levels and size bins equally.
%
% --- Parameter groups (from Adrian, June 2026) ---
%   Group 1 — fit first:
%     alpha (stickiness, 0.1 to 1.0)   <- this script
%     fp_alpha_cross
%     mining_dm
%   Group 2 — fit second:
%     fr_dim, zoo_p
%   Group 3 — fix:
%     fp_excess_density, settling law constants
%
% Steps:
%   1. Build real DepthProfile from keps_for_dave.mat
%   2. Build UVP observed phi(z): cruise-mean, mapped to model bins
%   3. Build daily surface forcing from UVP surface rows
%   4. Grid search: for each alpha, run daily-forced column, compute loss
%   5. Plot loss curve and best-fit depth profile

clear; close all; clc;

repo_root  = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir   = fullfile(repo_root, 'scripts', 'data');
addpath(genpath(fullfile(repo_root, 'src')));
addpath(data_dir);

% --- file paths ---
mat_path = fullfile(repo_root, 'data', 'NA', 'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(repo_root, 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% --- grid and config ---
col_grid = ColumnGrid(1000, 20);

base_cfg = SimulationConfig();
base_cfg.n_sections     = 30;
base_cfg.sinking_law    = 'kriest_8';
base_cfg.disagg_mode    = 'operator_split';
base_cfg.disagg_dmax_cm = 1.0;      % fallback; depth run uses eps(z)
base_cfg.enable_coag    = true;
base_cfg.enable_sinking = true;
base_cfg.enable_disagg  = true;
base_cfg.enable_zoo     = true;
base_cfg.enable_microbe = true;
base_cfg.enable_mining  = true;
base_cfg.microbe_r0     = 0.001;
base_cfg.surface_pp_mu  = 0.1;
base_cfg.r_to_rg        = 1.6;
base_cfg.fp_alpha_cross = 0.5;
base_cfg.zoo_Zc         = 0.307;
base_cfg.zoo_c          = 0.025;    % Stemmann 2004
base_cfg.zoo_Zf         = 0.063;
base_cfg.zoo_s          = 1.3e-5;   % Stemmann 2004
base_cfg.zoo_p          = 0.5;      % Stemmann 2004
base_cfg.zoo_ic         = 7;
base_cfg.mining_s       = 1.3e-5;
base_cfg.validate();

dt            = 0.25;   % day
steps_per_day = round(1 / dt);

% =============================================================
% Step 1: real depth profile
% =============================================================
fprintf('Loading keps profile...\n');
prof = load_keps(mat_path, col_grid.z_centers);

% =============================================================
% Step 2: UVP observed cruise-mean phi(z) -> model bins
% =============================================================
fprintf('Parsing UVP file...\n');
uvp = parse_uvp(uvp_file);

% map UVP phi to model size bins (nearest by diameter)
grid_d    = base_cfg.derive();
r_cm      = (0.75 / pi * grid_d.av_vol(:)).^(1/3);
d_model   = 2 * r_cm * 1e4;   % um
n_sec     = base_cfg.n_sections;

bin_map = zeros(1, numel(uvp.d_um));
for i = 1:numel(uvp.d_um)
    [~, bin_map(i)] = min(abs(d_model - uvp.d_um(i)));
end

% sum aggregate-sized UVP phi into model bins at each depth
n_ud = numel(uvp.depth_m);
uvp_phi_bins = zeros(n_ud, n_sec);
for i = 1:numel(uvp.d_um)
    if uvp.d_um(i) >= 2000
        continue;
    end
    k = bin_map(i);
    vals = uvp.phi(:, i);
    vals(isnan(vals)) = 0;
    uvp_phi_bins(:, k) = uvp_phi_bins(:, k) + vals;
end

% interpolate UVP phi to model z grid (total phi for loss)
uvp_phi_total = sum(uvp_phi_bins, 2);
Y_obs_z = interp1(uvp.depth_m, uvp_phi_total, col_grid.z_centers, 'pchip', 'extrap');
Y_obs_z = max(Y_obs_z, 0);

% =============================================================
% Step 3: daily surface forcing
% =============================================================
fprintf('Building daily surface forcing...\n');
daily = get_daily_surface_phi(uvp_file, base_cfg, col_grid);
n_days = daily.n_days;
fprintf('  %d days (%d with real UVP data)\n', n_days, sum(daily.has_data));

% =============================================================
% Step 4: grid search
% =============================================================
alpha_vals = [0.01 0.02 0.05 0.075 0.1 : 0.1 : 1.0];
n_alpha    = length(alpha_vals);
losses     = zeros(n_alpha, 1);

% pre-build simulation objects (same structure, only alpha changes)
fprintf('\n--- Grid search over alpha ---\n');
tic;

for ia = 1:n_alpha
    cfg_i       = copy(base_cfg);
    cfg_i.alpha = alpha_vals(ia);

    sim_i = ColumnSimulation(cfg_i, col_grid, prof);

    n_z = col_grid.n_z;
    Y   = zeros(n_z, n_sec);
    Yfp = zeros(n_z, n_sec);

    Y_sum_depth = zeros(n_z, 1);

    % daily-forced run with true Dirichlet BC at surface
    for i_day = 1:n_days
        for i_step = 1:steps_per_day
            Y(1, :) = daily.phi(i_day, :);
            [Y, Yfp] = sim_i.rhs.stepY(Y, dt, Yfp);
            Y(1, :) = daily.phi(i_day, :);
        end

        Y_sum_depth = Y_sum_depth + sum(Y + Yfp, 2);
    end

    % model cruise-mean total phi vs UVP cruise-mean total phi
    Y_total = Y_sum_depth ./ n_days;   % n_z x 1

    losses(ia) = loss_size_dist_1d(Y_total, Y_obs_z);
    fprintf('  alpha = %.2f   loss = %.4e\n', alpha_vals(ia), losses(ia));
end
elapsed = toc;

% =============================================================
% Step 5: report and figures
% =============================================================
[L_best, i_best] = min(losses);
alpha_best = alpha_vals(i_best);
fprintf('\nBest alpha = %.2f  (loss = %.4e)\n', alpha_best, L_best);
fprintf('Total time = %.1f min\n', elapsed / 60);

% --- run best-alpha model one more time for profile plot ---
cfg_b       = copy(base_cfg);
cfg_b.alpha = alpha_best;
sim_b = ColumnSimulation(cfg_b, col_grid, prof);
n_z = col_grid.n_z;
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);
Y_best_sum = zeros(n_z, 1);
for i_day = 1:n_days
    for i_step = 1:steps_per_day
        Y(1, :) = daily.phi(i_day, :);
        [Y, Yfp] = sim_b.rhs.stepY(Y, dt, Yfp);
        Y(1, :) = daily.phi(i_day, :);
    end
    Y_best_sum = Y_best_sum + sum(Y + Yfp, 2);
end
Y_best_total = Y_best_sum ./ n_days;

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Figure 1: loss vs alpha
f1 = figure;
semilogy(alpha_vals, losses, 'ko-', 'MarkerSize', 5, 'LineWidth', 1.2);
hold on;
semilogy(alpha_best, L_best, 'r*', 'MarkerSize', 10);
hold off;
xlabel('alpha');
ylabel('loss');
title(sprintf('grid search: best alpha = %.2f', alpha_best));
saveas(f1, fullfile(fig_dir, 'opt_alpha_grid_loss.png'));

% Figure 2: best-fit depth profile vs UVP
f2 = figure;
semilogy(Y_best_total, col_grid.z_centers, 'b-',  'DisplayName', sprintf('model (\\alpha=%.2f)', alpha_best));
hold on;
semilogy(Y_obs_z,      col_grid.z_centers, 'r--', 'DisplayName', 'UVP observed');
hold off;
set(gca, 'YDir', 'reverse');
xlabel('\phi_{total}  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('best-fit alpha: model vs UVP');
saveas(f2, fullfile(fig_dir, 'opt_alpha_best_profile.png'));

fprintf('\nFigures saved to docs/figures/\n');

% =============================================================
% Local loss function
% =============================================================
function L = loss_size_dist_1d(Y_model, Y_obs)
% Log-space squared difference between model and observed phi(z).
% Both inputs are n_z x 1 total phi vectors.
% Only include depths where both model and obs are > 0.
floor_val = 1e-15;
Y_m = max(Y_model(:), floor_val);
Y_o = max(Y_obs(:),   floor_val);
L = sum((log10(Y_m) - log10(Y_o)).^2);
end

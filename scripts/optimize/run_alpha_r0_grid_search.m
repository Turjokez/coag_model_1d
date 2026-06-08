% run_alpha_r0_grid_search.m
% 2-parameter grid search: alpha x microbe_r0, with spinup to steady state.
% Loss = sum of squared log10 differences, UVP <2000 um only (aggregates).

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir  = fullfile(repo_root, 'scripts', 'data');
addpath(genpath(fullfile(repo_root, 'src')));
addpath(data_dir);

mat_path = fullfile(repo_root, 'data', 'NA', 'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(repo_root, 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

col_grid = ColumnGrid(1000, 20);

base_cfg = SimulationConfig();
base_cfg.n_sections     = 30;
base_cfg.sinking_law    = 'kriest_8';
base_cfg.disagg_mode    = 'operator_split';
base_cfg.disagg_dmax_cm = 1.0;
base_cfg.enable_coag    = true;
base_cfg.enable_sinking = true;
base_cfg.enable_disagg  = true;
base_cfg.enable_zoo     = true;
base_cfg.enable_microbe = true;
base_cfg.enable_mining  = true;
base_cfg.surface_pp_mu  = 0.1;
base_cfg.r_to_rg        = 1.6;
base_cfg.fp_alpha_cross = 0.5;
base_cfg.zoo_Zc         = 0.307;
base_cfg.zoo_Zf         = 0.063;
base_cfg.zoo_ic         = 7;

dt            = 0.25;          % day
steps_per_day = round(1/dt);
max_cycles    = 20;
spinup_tol    = 0.01;

% load real profiles
prof = load_keps(mat_path, col_grid.z_centers);

% UVP observed phi(z), aggregates only (< 2000 um)
uvp = parse_uvp(uvp_file);

grid_d  = base_cfg.derive();
r_cm    = (0.75/pi * grid_d.av_vol(:)).^(1/3);
d_model = 2*r_cm*1e4;           % model bin diameters [um]
n_sec   = base_cfg.n_sections;

d_max_um = 2000;                 % exclude zooplankton-dominated bins
uvp_mask = uvp.d_um < d_max_um;

bin_map = zeros(1, numel(uvp.d_um));
for i = 1:numel(uvp.d_um)
    [~, bin_map(i)] = min(abs(d_model - uvp.d_um(i)));
end

n_ud = numel(uvp.depth_m);
uvp_phi_bins = zeros(n_ud, n_sec);
for i = 1:numel(uvp.d_um)
    if ~uvp_mask(i), continue; end
    k = bin_map(i);
    v = uvp.phi(:,i);  v(isnan(v)) = 0;
    uvp_phi_bins(:,k) = uvp_phi_bins(:,k) + v;
end

uvp_total = sum(uvp_phi_bins, 2);
Y_obs_full = max(interp1(uvp.depth_m, uvp_total, col_grid.z_centers, 'pchip', 'extrap'), 0);

% clip comparison to top 500 m only
comp_depth = 500;
z_mask = col_grid.z_centers <= comp_depth;
Y_obs = Y_obs_full(z_mask);

% daily surface forcing, large bins zeroed
daily = get_daily_surface_phi(uvp_file, base_cfg, col_grid);
n_days = daily.n_days;
daily.phi(:, d_model > d_max_um) = 0;

% grid
alpha_vals = [0.01, 0.05, 0.1, 0.2, 0.5, 1.0];
r0_vals    = [0, 0.0001, 0.0005, 0.001, 0.005, 0.01];
n_a  = numel(alpha_vals);
n_r0 = numel(r0_vals);

loss_grid   = zeros(n_a, n_r0);
cycles_grid = zeros(n_a, n_r0);

n_z = col_grid.n_z;
fprintf('%d x %d runs with spinup...\n', n_a, n_r0);
tic;

for ia = 1:n_a
    for ir = 1:n_r0
        cfg_i            = copy(base_cfg);
        cfg_i.alpha      = alpha_vals(ia);
        cfg_i.microbe_r0 = r0_vals(ir);

        sim = ColumnSimulation(cfg_i, col_grid, prof);
        Y   = zeros(n_z, n_sec);
        Yfp = zeros(n_z, n_sec);

        for ic = 1:max_cycles
            prev = mean(sum(Y + Yfp, 2));
            for id = 1:n_days
                for is = 1:steps_per_day
                    Y(1,:) = daily.phi(id,:);
                    [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
                    Y(1,:) = daily.phi(id,:);
                end
            end
            cur = mean(sum(Y + Yfp, 2));
            if prev > 0 && abs(cur-prev)/prev < spinup_tol && ic > 2, break; end
        end
        cycles_grid(ia, ir) = ic;
        loss_grid(ia, ir)   = log_loss(sum(Y+Yfp, 2), Y_obs, z_mask);
    end
    fprintf('  alpha=%.3f  cycles=%d-%d\n', alpha_vals(ia), ...
        min(cycles_grid(ia,:)), max(cycles_grid(ia,:)));
end
fprintf('done in %.1f min\n', toc/60);

[L_best, idx] = min(loss_grid(:));
[ia_best, ir_best] = ind2sub([n_a, n_r0], idx);
alpha_best = alpha_vals(ia_best);
r0_best    = r0_vals(ir_best);
fprintf('best: alpha=%.3f  r0=%.4f  loss=%.4e\n', alpha_best, r0_best, L_best);

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% loss surface
figure;
imagesc(log10(r0_vals + 1e-10), alpha_vals, loss_grid);
colorbar;
hold on;
plot(log10(r0_best + 1e-10), alpha_best, 'r+', 'MarkerSize', 12, 'LineWidth', 2);
hold off;
xlabel('log_{10}(microbe\_r0)');
ylabel('\alpha');
title('loss: \alpha x r_0');
saveas(gcf, fullfile(fig_dir, 'opt_2d_loss_surface.png'));

% best-fit profile: run to convergence
cfg_b            = copy(base_cfg);
cfg_b.alpha      = alpha_best;
cfg_b.microbe_r0 = r0_best;
sim_b = ColumnSimulation(cfg_b, col_grid, prof);
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);
for ic = 1:max_cycles
    prev = mean(sum(Y+Yfp, 2));
    for id = 1:n_days
        for is = 1:steps_per_day
            Y(1,:) = daily.phi(id,:);
            [Y, Yfp] = sim_b.rhs.stepY(Y, dt, Yfp);
            Y(1,:) = daily.phi(id,:);
        end
    end
    cur = mean(sum(Y+Yfp, 2));
    if prev > 0 && abs(cur-prev)/prev < spinup_tol && ic > 2, break; end
end

figure;
semilogy(sum(Y+Yfp,2), col_grid.z_centers, 'b-', ...
    'DisplayName', sprintf('model (\\alpha=%.2f, r_0=%.4f)', alpha_best, r0_best));
hold on;
semilogy(Y_obs_full, col_grid.z_centers, 'r--', 'DisplayName', 'UVP <2mm');
yline(comp_depth, 'k:', '500 m');
hold off;
set(gca, 'YDir', 'reverse');
xlabel('\phi  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('best fit vs UVP');
saveas(gcf, fullfile(fig_dir, 'opt_2d_best_profile.png'));

function L = log_loss(Y_model, Y_obs, z_mask)
% log-space loss over the comparison depth range only
Ym = sum(Y_model(z_mask, :), 2);
f = 1e-15;
L = sum((log10(max(Ym(:),f)) - log10(max(Y_obs(:),f))).^2);
end

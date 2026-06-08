
% run_spinup.m
% Repeat the 26-day UVP forcing cycle until column reaches steady state.
% Compares converged model phi(z) to UVP cruise-mean phi(z), <2000 um only.

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
data_dir  = fullfile(repo_root, 'scripts', 'data');
addpath(genpath(fullfile(repo_root, 'src')));
addpath(data_dir);

mat_path = fullfile(repo_root, 'data', 'NA', 'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(repo_root, 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

alpha_test = 0.1;

col_grid = ColumnGrid(1000, 20);

cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.enable_coag    = true;
cfg.enable_sinking = true;
cfg.enable_disagg  = true;
cfg.enable_zoo     = true;
cfg.enable_microbe = true;
cfg.enable_mining  = true;
cfg.microbe_r0     = 0.001;
cfg.surface_pp_mu  = 0.1;
cfg.r_to_rg        = 1.6;
cfg.fp_alpha_cross = 0.5;
cfg.zoo_Zc         = 0.307;
cfg.zoo_Zf         = 0.063;
cfg.zoo_ic         = 7;
cfg.alpha          = alpha_test;

dt            = 0.25;
steps_per_day = round(1/dt);
max_cycles    = 20;
tol           = 0.01;
d_max_um      = 2000;           % exclude zooplankton-dominated bins

prof  = load_keps(mat_path, col_grid.z_centers);
daily = get_daily_surface_phi(uvp_file, cfg, col_grid);
n_days = daily.n_days;

% UVP observed, aggregates only
uvp = parse_uvp(uvp_file);
grid_d  = cfg.derive();
r_cm    = (0.75/pi * grid_d.av_vol(:)).^(1/3);
d_model = 2*r_cm*1e4;
n_sec   = cfg.n_sections;

bin_map = zeros(1, numel(uvp.d_um));
for i = 1:numel(uvp.d_um)
    [~, bin_map(i)] = min(abs(d_model - uvp.d_um(i)));
end
n_ud = numel(uvp.depth_m);
uvp_phi_bins = zeros(n_ud, n_sec);
for i = 1:numel(uvp.d_um)
    if uvp.d_um(i) >= d_max_um, continue; end
    k = bin_map(i);
    v = uvp.phi(:,i);  v(isnan(v)) = 0;
    uvp_phi_bins(:,k) = uvp_phi_bins(:,k) + v;
end
uvp_total = sum(uvp_phi_bins, 2);
Y_obs = max(interp1(uvp.depth_m, uvp_total, col_grid.z_centers, 'pchip', 'extrap'), 0);

comp_depth = 500;                   % compare only top 500 m

% also fix the daily forcing filter correctly
daily.phi(:, d_model > d_max_um) = 0;

% spinup
sim = ColumnSimulation(cfg, col_grid, prof);
n_z = col_grid.n_z;
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);

profiles   = zeros(max_cycles, n_z);
conv_trace = zeros(max_cycles, 1);

fprintf('spinup: alpha=%.3f\n', alpha_test);
tic;
for ic = 1:max_cycles
    prev = mean(sum(Y+Yfp, 2));
    for id = 1:n_days
        for is = 1:steps_per_day
            Y(1,:) = daily.phi(id,:);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(1,:) = daily.phi(id,:);
        end
    end
    cur = mean(sum(Y+Yfp, 2));
    if prev > 0
        conv_trace(ic) = abs(cur-prev)/prev;
    else
        conv_trace(ic) = 1;
    end
    profiles(ic,:) = sum(Y+Yfp, 2)';
    fprintf('  cycle %2d: depth-mean=%.3e  rel change=%.3f\n', ic, cur, conv_trace(ic));
    if conv_trace(ic) < tol && ic > 2
        fprintf('converged at cycle %d\n', ic);
        break;
    end
end
n_done = ic;
fprintf('time: %.1f min\n', toc/60);

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% depth profiles by cycle (show full column, mark comparison limit)
figure;
cm = parula(n_done);
for ic = 1:n_done
    semilogy(profiles(ic,:)', col_grid.z_centers, '-', 'Color', cm(ic,:));
    hold on;
end
semilogy(Y_obs, col_grid.z_centers, 'r--', 'LineWidth', 1.5);
yline(comp_depth, 'k:', '500 m');
set(gca, 'YDir', 'reverse');
xlabel('\phi  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
c = colorbar;  c.Label.String = 'cycle';
c.Ticks = [0 1];  c.TickLabels = {'1', num2str(n_done)};
title(sprintf('spinup \\alpha=%.2f (%d cycles)', alpha_test, n_done));
saveas(gcf, fullfile(fig_dir, 'spinup_profiles.png'));

% convergence
figure;
semilogy(1:n_done, conv_trace(1:n_done), 'ko-', 'MarkerSize', 5);
hold on;
yline(tol, 'r--');
xlabel('cycle');
ylabel('rel change');
title('spinup convergence');
saveas(gcf, fullfile(fig_dir, 'spinup_convergence.png'));

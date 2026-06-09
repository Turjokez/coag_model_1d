% run_data_bc_compare.m
%
% Compare two surface boundary choices:
%   1. reset before each substep only
%   2. reset before and after each substep
%
% This checks if the subsurface peak is caused by the surface reset rule.

clear; close all; clc;

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

cfg = SimulationConfig();
cfg.n_sections      = 30;
cfg.sinking_law     = 'kriest_8';
cfg.disagg_mode     = 'operator_split';
cfg.disagg_dmax_cm  = 1.0;      % fallback; depth run uses eps(z)
cfg.enable_coag     = true;
cfg.enable_disagg   = true;
cfg.enable_zoo      = true;
cfg.enable_microbe  = true;
cfg.enable_mining   = true;
cfg.microbe_r0      = 0.001;
cfg.r_to_rg         = 1.6;
cfg.zoo_c           = 0.025;    % Stemmann 2004
cfg.zoo_s           = 1.3e-5;   % Stemmann 2004
cfg.zoo_p           = 0.5;      % Stemmann 2004
cfg.zoo_ic          = 7;        % bin 8, about 115 um
cfg.mining_s        = 1.3e-5;
cfg.fp_alpha_cross  = 0.5;
cfg.validate();

dt = 0.25;                  % day
steps_per_day = round(1 / dt);

col_grid = ColumnGrid(1000, 20);
prof  = load_keps(mat_path, col_grid.z_centers);
daily = get_daily_surface_phi(uvp_file, cfg, col_grid);
uvp   = parse_uvp(uvp_file);

% UVP cruise-mean aggregate phi on model depths
mask_agg = uvp.d_um < 2000;
uvp_phi = uvp.phi(:, mask_agg);
uvp_phi(isnan(uvp_phi)) = 0;
uvp_total = sum(uvp_phi, 2);
uvp_model = interp1(uvp.depth_m, uvp_total, col_grid.z_centers, 'pchip', 'extrap');
uvp_model = max(0, uvp_model);

fprintf('Running boundary condition comparison...\n');
[phi_pre, surf_pre]   = run_case(false, cfg, col_grid, prof, daily, dt, steps_per_day);
[phi_both, surf_both] = run_case(true,  cfg, col_grid, prof, daily, dt, steps_per_day);

[~, iz200] = min(abs(col_grid.z_centers - 200));
fprintf('\n--- Summary ---\n');
fprintf('UVP surface mean:        %.3e\n', mean(sum(daily.phi, 2)));
fprintf('Pre-only surface mean:   %.3e\n', mean(surf_pre));
fprintf('Both-reset surface mean: %.3e\n', mean(surf_both));
fprintf('UVP %.0f m phi:          %.3e\n', col_grid.z_centers(iz200), uvp_model(iz200));
fprintf('Pre-only %.0f m phi:     %.3e\n', col_grid.z_centers(iz200), phi_pre(iz200));
fprintf('Both-reset %.0f m phi:   %.3e\n', col_grid.z_centers(iz200), phi_both(iz200));

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

figure;
semilogy(phi_pre,   col_grid.z_centers, 'b-',  'DisplayName', 'pre only');
hold on;
semilogy(phi_both,  col_grid.z_centers, 'k-',  'DisplayName', 'pre + post');
semilogy(uvp_model, col_grid.z_centers, 'r--', 'DisplayName', 'UVP <2000 um');
hold off;
set(gca, 'YDir', 'reverse');
xlabel('\phi_{total}  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('surface reset test');
saveas(gcf, fullfile(fig_dir, 'data_bc_compare_depth.png'));

figure;
plot(daily.day_num, sum(daily.phi, 2), 'r-o', 'MarkerSize', 3, ...
    'DisplayName', 'UVP forcing');
hold on;
plot(daily.day_num, surf_pre,  'b-', 'DisplayName', 'pre only');
plot(daily.day_num, surf_both, 'k-', 'DisplayName', 'pre + post');
hold off;
xlabel('day');
ylabel('\phi  [cm^3 cm^{-3}]');
legend('location', 'best');
title('surface boundary check');
saveas(gcf, fullfile(fig_dir, 'data_bc_compare_surface.png'));

fprintf('\nSaved figures:\n');
fprintf('  docs/figures/data_bc_compare_depth.png\n');
fprintf('  docs/figures/data_bc_compare_surface.png\n');

function [phi_mean, surf_total] = run_case(reset_after, cfg, col_grid, prof, daily, dt, steps_per_day)
% Run one boundary condition case.
n_z   = col_grid.n_z;
n_sec = cfg.n_sections;
n_days = daily.n_days;

sim = ColumnSimulation(cfg, col_grid, prof);

Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);
Y(1, :) = daily.phi(1, :);

phi_day = zeros(n_days, n_z);
surf_total = zeros(n_days, 1);

for i_day = 1:n_days
    for i_step = 1:steps_per_day
        Y(1, :) = daily.phi(i_day, :);
        [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
        if reset_after
            Y(1, :) = daily.phi(i_day, :);
        end
    end

    total = Y + Yfp;
    phi_day(i_day, :) = sum(total, 2)';
    surf_total(i_day) = sum(total(1, :));
end

phi_mean = mean(phi_day, 1)';
end

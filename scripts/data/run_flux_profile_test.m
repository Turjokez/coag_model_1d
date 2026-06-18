% run_flux_profile_test.m
%
% Test A: Flux profile and Transfer Efficiency.
%
% Model was designed to compute sinking flux, not standing stock.
% This script asks: even if model/UVP standing stock ratio is off,
% does the MODEL FLUX profile have the right shape?
%
% Flux at depth z:  F(z) = sum_bins( w_bin * Y(z, bin) )   [m/day * BV_units]
%
% We compare:
%   1. Model flux profile F(z) vs Martin power law: F(z) = F(z0) * (z/z0)^(-b)
%      Martin b ~ 0.85-1.0 for open ocean
%   2. Transfer efficiency: TE = F(500m) / F(100m)
%   3. Flux attenuation slope on log-log axes
%
% Uses flux BC (more physical than Dirichlet for flux interpretation).
% Also runs Dirichlet for comparison.
%
% No external flux data needed. This is a self-consistency check.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Setup
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

k_bc      = 2;
dz        = col_grid.dz;
z_centers = col_grid.z_centers;   % all 20 layers
n_z       = col_grid.n_z;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
k_plot_bc = 2:10;   % for BC loading only
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot_bc);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% sinking speeds: Kriest_8
d_cm  = bc.d_model_um * 1e-4;
w_bin = (66 * d_cm .^ 0.62)';   % 1 x n_sec [m/day]

% reference depth indices
[~, k100]  = min(abs(z_centers - 100));   % ~100 m
[~, k500]  = min(abs(z_centers - 500));   % ~500 m
[~, k1000] = min(abs(z_centers - 975));   % bottom

[~, ia, ib] = intersect(bc.dates, uvpd.dates);

% ---------------------------------------------------------------
% 2. Two cases: Dirichlet and Flux BC
% ---------------------------------------------------------------
cases = { ...
    struct('bc_type', 'dirichlet', 'label', 'Dirichlet BC', 'color', 'k', 'ls', '-'), ...
    struct('bc_type', 'flux',      'label', 'Flux BC',      'color', 'b', 'ls', '--'), ...
};

flux_profiles = zeros(n_z, 2);   % mean flux at each layer

for ic = 1:2
    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    use_flux = strcmp(cases{ic}.bc_type, 'flux');

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
            for i_step = 1:steps_per_day
                if use_flux
                    Y(k_bc, :) = Y(k_bc, :) + flux_src;
                    [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                else
                    Y(k_bc, :) = phi_bc_daily(i_day, :);
                    [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                    Y(k_bc, :) = phi_bc_daily(i_day, :);
                end
            end
        end
        phi_after  = mean(sum(Y + Yfp, 2));
        rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
        if rel_change < spinup_tol
            fprintf('%s: converged at cycle %d\n', cases{ic}.label, icyc);
            break;
        end
    end

    % final run: accumulate flux on cast days
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    flux_acc = zeros(n_z, 1);
    n_cast = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            if use_flux
                Y(k_bc, :) = Y(k_bc, :) + flux_src;
                [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            else
                Y(k_bc, :) = phi_bc_daily(i_day, :);
                [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                Y(k_bc, :) = phi_bc_daily(i_day, :);
            end
        end
        if any(bc.dates(i_day) == uvpd.dates)
            % flux at each layer = w * (Y + Yfp)
            Ytot = Y + Yfp;
            for k = 1:n_z
                flux_acc(k) = flux_acc(k) + sum(w_bin .* Ytot(k, :));
            end
            n_cast = n_cast + 1;
        end
    end

    flux_profiles(:, ic) = flux_acc / max(n_cast, 1);
end

% ---------------------------------------------------------------
% 3. Martin curve reference
% ---------------------------------------------------------------
% Martin (1987): F(z) = F(z0) * (z/z0)^(-b), b = 0.858
% We normalize to model flux at k100 (Dirichlet case)
z0   = z_centers(k100);
b_martin = 0.858;
F_ref = flux_profiles(k100, 1);   % normalize to Dirichlet at 100m
z_martin = z_centers(k100:end);
F_martin = F_ref * (z_martin / z0) .^ (-b_martin);

% ---------------------------------------------------------------
% 4. Print transfer efficiency
% ---------------------------------------------------------------
fprintf('\n--- Transfer Efficiency: F(z) / F(100m) ---\n');
fprintf('%-10s  %-16s  %-10s\n', 'Depth', cases{1}.label, cases{2}.label);
check_depths = [200 300 500 750 975];
for zd = check_depths
    [~, kz] = min(abs(z_centers - zd));
    te1 = flux_profiles(kz, 1) / max(flux_profiles(k100, 1), 1e-30);
    te2 = flux_profiles(kz, 2) / max(flux_profiles(k100, 2), 1e-30);
    fprintf('%5.0f m     %5.3f             %5.3f\n', z_centers(kz), te1, te2);
end

% Martin b exponent (log-log slope between 100m and 500m)
for ic = 1:2
    F_top  = flux_profiles(k100, ic);
    F_bot  = flux_profiles(k500, ic);
    z_top  = z_centers(k100);
    z_bot  = z_centers(k500);
    if F_top > 0 && F_bot > 0
        b_model = -log(F_bot / F_top) / log(z_bot / z_top);
        fprintf('%s: Martin b (100-500m) = %.3f\n', cases{ic}.label, b_model);
    end
end

% ---------------------------------------------------------------
% 5. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 14], 'Color', 'white');
hold on;

% Martin reference
semilogy(F_martin, z_martin, 'r:', 'LineWidth', 1.0, 'DisplayName', ...
    sprintf('Martin b=%.2f', b_martin));

for ic = 1:2
    semilogy(flux_profiles(:, ic), z_centers, ...
        [cases{ic}.color cases{ic}.ls 'o'], ...
        'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end

set(gca, 'YDir', 'reverse', 'YLim', [50 1000], 'YScale', 'linear');
xlabel('Flux  (w \times BV)');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Flux profile vs Martin curve', 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'flux_profile_test.png'));
fprintf('\nSaved flux_profile_test.png\n');

% ---------------------------------------------------------------
function cfg = cfg_best()
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.enable_zoo     = true;
cfg.enable_microbe = false;
cfg.enable_mining  = true;
cfg.alpha          = 0.10;
cfg.microbe_r0     = 0.0;
cfg.surface_pp_mu  = 0.0;
cfg.r_to_rg        = 1.6;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.mining_s       = 1.3e-5;
cfg.fp_alpha_cross = 0.5;
end

% run_number_biovolume_diagnostic.m
%
% Adrian's diagnostic (June 12 2026 meeting):
% Calculate total particle number N(z) and biovolume BV(z) for
% model and UVP. Compare to distinguish two hypotheses:
%
%   DVM hypothesis:         UVP BV > model BV at depth (new mass)
%   Fragmentation hypothesis: UVP BV ~ model BV, but UVP N > model N
%
% Steps:
%   1. Run model (best config, 100m BC), save cast-day snapshots.
%   2. Load UVP data (100-2000 um).
%   3. Compute N and BV at each depth for each cast day.
%   4. Plot mean depth profiles: model vs UVP.

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
% 1. Config + grid (best config)
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

k_bc   = 2;
k_plot = 2:10;   % z = 75, 125, ..., 475 m
z_mod  = col_grid.z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.disagg_dmax_A  = 9.39e-6 * 5;   % Parker x5
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
cfg.validate();

n_sec  = cfg.n_sections;
grid_c = cfg.derive();
av_vol = grid_c.av_vol(:);   % cm^3 per particle per bin [n_sec x 1]

% model bin diameters [um] — to filter same range as UVP (100-2000 um)
r_mod_cm   = (0.75 / pi * av_vol).^(1/3);
d_mod_um   = 2 * r_mod_cm * 1e4;
mask_mod   = d_mod_um >= 100 & d_mod_um < 2000;   % match UVP range

% ---------------------------------------------------------------
% 2. BC + UVP
% ---------------------------------------------------------------
bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 3:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% UVP bins 100-2000 um only
mask_uvp   = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_um   = uvpd.d_um(mask_uvp);

% average particle volume per UVP bin [cm^3], assume sphere
r_uvp_cm   = (d_uvp_um / 2) / 1e4;   % um -> cm
av_vol_uvp = (4/3) * pi * r_uvp_cm.^3;   % cm^3

% UVP depth rows for our plot range
mask_z_uvp = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
z_uvp      = uvpd.depth_m(mask_z_uvp);

% ---------------------------------------------------------------
% 3. Spinup
% ---------------------------------------------------------------
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(col_grid.n_z, n_sec);
Yfp = zeros(col_grid.n_z, n_sec);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    if rel_change < spinup_tol
        fprintf('Spinup converged at cycle %d\n', icyc);
        break;
    end
end

% ---------------------------------------------------------------
% 4. Final run: save cast-day snapshots
% ---------------------------------------------------------------
Y   = zeros(col_grid.n_z, n_sec);
Yfp = zeros(col_grid.n_z, n_sec);
Y_daily = zeros(numel(k_plot), n_sec, n_days);

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    for i_step = 1:steps_per_day
        Y(k_bc, :) = phi_bc_daily(i_day, :);
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        Y(k_bc, :) = phi_bc_daily(i_day, :);
    end
    Ytot = Y + Yfp;
    Y_daily(:, :, i_day) = Ytot(k_plot, :);
end
fprintf('Model run complete\n');

% ---------------------------------------------------------------
% 5. Match cast days (model day index <-> UVP date index)
% ---------------------------------------------------------------
[~, ia, ib] = intersect(bc.dates, uvpd.dates);
n_cast = numel(ia);
fprintf('Cast days: %d\n', n_cast);

% ---------------------------------------------------------------
% 6. Compute N and BV at each depth, averaged over cast days
% ---------------------------------------------------------------
% Model: BV [cm^3/cm^3], N [particles/cm^3]
BV_mod_all = zeros(numel(k_plot), n_cast);
N_mod_all  = zeros(numel(k_plot), n_cast);

% UVP: BV [cm^3/cm^3], N [particles/cm^3]
BV_uvp_all = zeros(numel(z_uvp), n_cast);
N_uvp_all  = zeros(numel(z_uvp), n_cast);

for m = 1:n_cast
    id_mod = ia(m);
    id_uvp = ib(m);

    % model — BV uses all bins; N uses only 100-2000 um bins (same as UVP)
    phi_m = squeeze(Y_daily(:, :, id_mod));   % [n_z x n_sec]
    BV_mod_all(:, m) = sum(phi_m, 2);
    N_mod_all(:, m)  = sum(bsxfun(@rdivide, phi_m(:, mask_mod), av_vol(mask_mod)'), 2);

    % UVP
    phi_u = squeeze(uvpd.phi(id_uvp, mask_z_uvp, mask_uvp));   % [n_uvp_z x n_uvp_bins]
    if size(phi_u, 2) ~= numel(av_vol_uvp)
        phi_u = phi_u';
    end
    phi_u(isnan(phi_u)) = 0;
    BV_uvp_all(:, m) = sum(phi_u, 2);
    N_uvp_all(:, m)  = sum(bsxfun(@rdivide, phi_u, av_vol_uvp(:)'), 2);
end

% time-mean profiles
BV_mod = mean(BV_mod_all, 2);   % [n_z x 1]
N_mod  = mean(N_mod_all,  2);
BV_uvp = mean(BV_uvp_all, 2);   % [n_uvp_z x 1]
N_uvp  = mean(N_uvp_all,  2);

% ---------------------------------------------------------------
% 7. Plot: 2 panels — BV and N vs depth
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [1 1 16 10], 'Color', 'white');

% Panel a: Biovolume
ax1 = subplot(1, 2, 1);
plot(BV_uvp * 1e6, z_uvp, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
plot(BV_mod * 1e6, z_mod, 'r-s', 'MarkerSize', 3, 'LineWidth', 1.2);
set(ax1, 'YDir', 'reverse', 'YLim', [60 510], 'XScale', 'log');
xlabel('Biovolume (ppmV)', 'FontSize', 8);
ylabel('Depth (m)', 'FontSize', 8);
legend('UVP', 'Model', 'Location', 'southeast', 'FontSize', 7);
title('a) Total biovolume', 'FontSize', 8, 'FontWeight', 'normal');

% Panel b: Particle number
ax2 = subplot(1, 2, 2);
plot(N_uvp, z_uvp, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
plot(N_mod, z_mod, 'r-s', 'MarkerSize', 3, 'LineWidth', 1.2);
set(ax2, 'YDir', 'reverse', 'YLim', [60 510], 'XScale', 'log');
xlabel('Particle number (cm^{-3})', 'FontSize', 8);
legend('UVP', 'Model', 'Location', 'southeast', 'FontSize', 7);
title('b) Total particle number', 'FontSize', 8, 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'number_biovolume_diagnostic.png'));
fprintf('Saved number_biovolume_diagnostic.png\n');

% ---------------------------------------------------------------
% 8. Print ratios at each model depth
% ---------------------------------------------------------------
fprintf('\nDepth  BV_uvp/BV_mod  N_uvp/N_mod\n');
for k = 1:numel(z_mod)
    % interpolate UVP to model depth
    BV_u_k = interp1(z_uvp, BV_uvp, z_mod(k), 'linear', NaN);
    N_u_k  = interp1(z_uvp, N_uvp,  z_mod(k), 'linear', NaN);
    fprintf('%5.0f m   %6.2f          %6.2f\n', ...
        z_mod(k), BV_u_k / BV_mod(k), N_u_k / N_mod(k));
end

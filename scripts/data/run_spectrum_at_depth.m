% run_spectrum_at_depth.m
%
% Compare model vs UVP size spectrum at a specific target depth.
% Shows which size bins within 100-2000 um are most deficit in the model.
%
% Target depth: 375 m (center of deep residual problem).
% Config: alpha=0.10, Da x5, r0=0, 100m BC (best config).

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

target_z_m = 375;   % depth to compare [m]

% ---------------------------------------------------------------
% 1. Setup
% ---------------------------------------------------------------
cfg      = cfg_base();
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;
k_bc          = 2;

% find model layer closest to target depth
[~, k_target] = min(abs(col_grid.z_centers - target_z_m));
fprintf('Target depth: %d m -> model layer %d (z=%.0f m)\n', ...
    target_z_m, k_target, col_grid.z_centers(k_target));

bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 3:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% model bin diameters [um]
grid_c = cfg.derive();
av_vol = grid_c.av_vol(:);
r_cm   = (0.75 / pi * av_vol).^(1/3);
d_mod_um = 2 * r_cm * 1e4;
mask_mod = d_mod_um >= 100 & d_mod_um < 2000;

% UVP bin diameters and bins at target depth
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_um = uvpd.d_um(mask_uvp);
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
z_uvp    = uvpd.depth_m(mask_z);
[~, iz_target] = min(abs(z_uvp - target_z_m));

% ---------------------------------------------------------------
% 2. Spinup
% ---------------------------------------------------------------
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);

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
% 3. Final run — average over cast days
% ---------------------------------------------------------------
Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);
phi_mod_sum = zeros(1, cfg.n_sections);
phi_uvp_sum = zeros(1, sum(mask_uvp));
n_cast = 0;

[~, ia, ib] = intersect(bc.dates, uvpd.dates);

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    for i_step = 1:steps_per_day
        Y(k_bc, :) = phi_bc_daily(i_day, :);
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        Y(k_bc, :) = phi_bc_daily(i_day, :);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Ytot        = Y + Yfp;
        phi_mod_sum = phi_mod_sum + Ytot(k_target, :);

        % match this cast day to UVP
        [~, loc] = ismember(bc.dates(i_day), uvpd.dates);
        if loc > 0
            phi_u = squeeze(uvpd.phi(loc, mask_z, mask_uvp));
            if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
            phi_uvp_sum = phi_uvp_sum + phi_u(iz_target, :);
        end
        n_cast = n_cast + 1;
    end
end
phi_mod_mean = phi_mod_sum / max(n_cast, 1);
phi_uvp_mean = phi_uvp_sum / max(n_cast, 1);

% only UVP-range bins for model
phi_mod_mid = phi_mod_mean(mask_mod);
d_mod_mid   = d_mod_um(mask_mod);

% ---------------------------------------------------------------
% 4. Print bin-by-bin comparison
% ---------------------------------------------------------------
fprintf('\nModel vs UVP spectrum at ~%d m (cast-day mean)\n', ...
    col_grid.z_centers(k_target));

% interpolate UVP to model bin centres for ratio
phi_uvp_interp = interp1(log(d_uvp_um), phi_uvp_mean(:), log(d_mod_mid), 'linear', NaN);

fprintf('d_mod (um)   phi_mod      phi_uvp(interp)   ratio\n');
for i = 1:numel(d_mod_mid)
    fprintf('  %6.0f    %.3e      %.3e          %.2f\n', ...
        d_mod_mid(i), phi_mod_mid(i), phi_uvp_interp(i), ...
        phi_mod_mid(i) / max(phi_uvp_interp(i), 1e-30));
end

% ---------------------------------------------------------------
% 5. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 8], 'Color', 'white');
loglog(d_uvp_um, phi_uvp_mean, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
loglog(d_mod_mid, phi_mod_mid, 'r-s', 'MarkerSize', 4, 'LineWidth', 1.2);
xlabel('Diameter (\mum)');
ylabel('\phi (ppmV)');
legend('UVP', 'Model', 'Location', 'northeast', 'FontSize', 7);
title(sprintf('Size spectrum at %d m', col_grid.z_centers(k_target)), ...
    'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'spectrum_at_depth.png'));
fprintf('\nSaved spectrum_at_depth.png\n');

% ---------------------------------------------------------------
function cfg = cfg_base()
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

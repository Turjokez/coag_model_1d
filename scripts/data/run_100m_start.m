% run_100m_start.m
%
% Start the model at 100 m using UVP data there as the top BC.
% Predict what happens from 100 to 500 m and compare to UVP casts.
%
% Adrian (June 11): "below 100 m the main driver is what sinks from above.
% It is a cleaner test. Start there first."
%
% Strategy:
%   - Column: 1000 m, 20 layers, dz = 50 m (same grid as always).
%   - Top BC: layer 2 (z = 75 m, closest to 100 m) is reset to UVP at 100 m
%     each time step instead of surface UVP.
%   - Layer 1 is left free (model fills it; we ignore it in comparison).
%   - Comparison: layers 3-10 (z = 125-525 m) vs UVP cast data.
%
% Output: two figures
%   (1) phi profile: model vs UVP, 100-500 m
%   (2) size spectrum panels at selected depths (same as report figure)

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% ---------------------------------------------------------------
% 1. Config  (same as mass budget run)
% ---------------------------------------------------------------
cfg = SimulationConfig();
cfg.n_sections       = 30;
cfg.sinking_law      = 'kriest_8';
cfg.disagg_mode      = 'operator_split';
cfg.disagg_dmax_cm   = 1.0;
cfg.disagg_dmax_A    = 9.39e-6 * 5;
cfg.enable_coag      = true;
cfg.enable_disagg    = true;
cfg.enable_zoo       = true;
cfg.enable_microbe   = true;
cfg.enable_mining    = true;
cfg.alpha            = 0.5;
cfg.microbe_r0       = 0.03;
cfg.microbe_use_temp = true;
cfg.microbe_tref_C   = 20;
cfg.surface_pp_mu    = 0.0;   % no surface production below 100 m
cfg.r_to_rg          = 1.6;
cfg.zoo_c            = 0.025;
cfg.zoo_s            = 1.3e-5;
cfg.zoo_p            = 0.5;
cfg.zoo_ic           = 7;
cfg.mining_s         = 1.3e-5;
cfg.fp_alpha_cross   = 0.5;
cfg.validate();

col_grid  = ColumnGrid(1000, 20);      % dz = 50 m, centers at 25, 75, 125 ... 975 m
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);

% layer index for 100 m BC (layer 2, center = 75 m)
k_bc = 2;
fprintf('BC layer: %d  (z = %.0f m)\n', k_bc, col_grid.z_centers(k_bc));

% layers used for comparison (125-525 m)
k_compare = 3:10;
z_compare = col_grid.z_centers(k_compare);
fprintf('Compare layers: %d to %d  (z = %.0f to %.0f m)\n', ...
    k_compare(1), k_compare(end), z_compare(1), z_compare(end));

% ---------------------------------------------------------------
% 2. Build 100 m BC with power-law fill below 100 um
% ---------------------------------------------------------------
bc = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_compare);
uvpd           = bc.uvpd;
iz_100         = bc.iz_bc;
d_model_um     = bc.d_model_um;
mask_uvp_model = bc.mask_uvp_model;
phi_bc_daily   = bc.phi_bc_daily;
n_days         = bc.n_days;
n_sec          = cfg.n_sections;
id_model_best  = bc.id_model_best;
id_uvp_best    = bc.id_uvp_best;
phi_uvp_spec   = bc.phi_uvp_spec;
d_uvp_ok       = bc.d_uvp_ok;
dw_uvp_ok      = bc.dw_uvp_ok;
best_date      = bc.best_date;

fprintf('UVP depth used for BC: %.1f m (index %d)\n', bc.bc_depth_m, iz_100);
fprintf('BC phi stats (non-zero rows): %d of %d days have UVP cast at 100m\n', ...
    sum(any(phi_bc_daily > 0, 2)), n_days);
fprintf('Best cast day: %d (model index %d, UVP index %d)\n', ...
    bc.best_date, id_model_best, id_uvp_best);

% ---------------------------------------------------------------
% 4. Spinup: BC at layer k_bc = 2 (not surface)
% ---------------------------------------------------------------
dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(col_grid.n_z, n_sec);
Yfp = zeros(col_grid.n_z, n_sec);

fprintf('\nRunning spinup...\n');
for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);   % BC at 100 m
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    fprintf('  cycle %d: rel_change = %.4f\n', icyc, rel_change);
    if rel_change < spinup_tol
        fprintf('Converged at cycle %d\n', icyc);
        break;
    end
end

% second pass: capture snapshot on best cast day
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

% ---------------------------------------------------------------
% 5. UVP data at comparison depths on best cast day
% ---------------------------------------------------------------
% total phi per depth: model (UVP-range bins) vs UVP
phi_mod_cmp = sum(Y_snap(k_compare, mask_uvp_model), 2);   % n_cmp x 1
phi_uvp_tot = bc.phi_uvp_cmp;                              % n_cmp x 1
n_cmp = numel(k_compare);

% ---------------------------------------------------------------
% 6. Figure 1: total phi profile, 100-500 m
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

figure('Units','centimeters','Position',[2 2 9 13]);
hold on;
plot(phi_mod_cmp * 1e6, z_compare, 'b-', 'LineWidth', 1.5, ...
    'DisplayName', 'model (UVP range)');
plot(phi_uvp_tot * 1e6, z_compare, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'UVP');
% mark BC depth
yline(uvpd.depth_m(iz_100), 'r:', 'BC', 'LabelHorizontalAlignment','left');
set(gca, 'YDir', 'reverse');
xlabel('\phi [ppmV]');
ylabel('Depth (m)');
ylim([75 550]);
legend('location', 'southeast');
title(sprintf('100m start: %d', best_date));

saveas(gcf, fullfile(fig_dir, '100m_phi_profile.png'));
fprintf('Saved 100m_phi_profile.png\n');

% ---------------------------------------------------------------
% 7. Figure 2: size spectrum panels at selected depths
% ---------------------------------------------------------------
% show 3 panels: 125 m, 225 m, 425 m (layers 3, 5, 8 in k_compare)
panel_idx = [1, 3, 6];   % within k_compare  (z ~ 125, 225, 425 m)
n_panels  = numel(panel_idx);

% model spectrum: S = phi / delta_d [ppmV / mm]
d_model_edges = zeros(1, n_sec + 1);
d_model_edges(1)       = d_model_um(1)^2 / d_model_um(2);
d_model_edges(n_sec+1) = d_model_um(n_sec)^2 / d_model_um(n_sec-1);
for k = 2:n_sec
    d_model_edges(k) = sqrt(d_model_um(k-1) * d_model_um(k));
end
dw_model_mm = diff(d_model_edges) / 1000;  % um -> mm
d_model_mm  = d_model_um / 1000;

% UVP spectrum: S = phi / delta_d_mm [ppmV / mm], only UVP-range bins
dw_uvp_mm  = dw_uvp_ok / 1000;
d_uvp_mm   = d_uvp_ok  / 1000;

figure('Units','centimeters','Position',[2 2 18 6]);
for ip = 1:n_panels
    ki  = panel_idx(ip);   % index within k_compare
    k   = k_compare(ki);   % absolute layer
    z_m = z_compare(ki);

    S_mod = Y_snap(k, :) * 1e6 ./ dw_model_mm(:)';
    S_uvp = phi_uvp_spec(ki, :) * 1e6 ./ dw_uvp_mm(:)';

    subplot(1, n_panels, ip);
    hold on;
    plot(d_model_mm, S_mod, 'b-', 'LineWidth', 1.5);
    plot(d_uvp_mm,   S_uvp, 'k.', 'MarkerSize', 6);
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([0.05 5]);
    xlabel('d (mm)');
    if ip == 1
        ylabel('S [ppmV mm^{-1}]');
        legend('model', 'UVP', 'location', 'northeast');
    end
    title(sprintf('z = %d m', round(z_m)));
end

saveas(gcf, fullfile(fig_dir, '100m_spectrum_panels.png'));
fprintf('Saved 100m_spectrum_panels.png\n');

% ---------------------------------------------------------------
% 8. Print comparison table
% ---------------------------------------------------------------
fprintf('\n--- Model vs UVP: 100m start (best cast %d) ---\n', ...
    best_date);
fprintf('%-8s  %-12s  %-12s  %-8s\n', ...
    'Depth(m)', 'Model', 'UVP', 'Ratio');
fprintf('%-8s  %-12s  %-12s  %-8s\n', '', '[ppmV]', '[ppmV]', 'Mod/UVP');
fprintf('%s\n', repmat('-', 1, 45));
for i = 1:n_cmp
    rat = phi_mod_cmp(i) / max(phi_uvp_tot(i), 1e-30);
    fprintf('%8.1f  %12.3f  %12.3f  %8.2f\n', ...
        z_compare(i), phi_mod_cmp(i)*1e6, phi_uvp_tot(i)*1e6, rat);
end

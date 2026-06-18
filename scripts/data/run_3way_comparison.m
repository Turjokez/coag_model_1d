% run_3way_comparison.m
%
% Three-way model validation:
%   1. UVP standing stock: model BV / UVP BV at each depth
%   2. Flux profile: model flux vs Martin curve (b = 0.858)
%   3. Sediment trap: model aggregate flux vs Durkin gel trap
%
% All three use flux BC (more physical than Dirichlet).
% Output: 2-panel figure + summary table.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
trap_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'sediment_trap_durkin', 'raw', ...
    'cb6a494508_EXPORTS_EXPORTSNA_JC214_classified_geltrap_particlefluxes.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Setup
% ---------------------------------------------------------------
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);
z_centers = col_grid.z_centers;

k_bc   = 2;
dz     = col_grid.dz;
n_z    = col_grid.n_z;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
k_plot = 2:10;
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;
d_model_um   = bc.d_model_um;

% sinking speeds and particle volumes
d_cm  = d_model_um * 1e-4;
d_m   = d_model_um * 1e-6;
w_bin = (66 * d_cm .^ 0.62)';     % 1 x n_sec [m/day]
V_bin = ((pi/6) * d_m .^ 3)';     % 1 x n_sec [m3/particle]

mask_uvp = d_model_um >= 100 & d_model_um < 2000;

% ---------------------------------------------------------------
% 2. UVP reference (standing stock)
% ---------------------------------------------------------------
mask_uvp_raw = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z_uvp   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
z_mod        = z_centers(k_plot);
[~, ia, ib]  = intersect(bc.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    phi_u  = squeeze(uvpd.phi(id_uvp, mask_z_uvp, mask_uvp_raw));
    if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z_uvp) - z_mod(ki)));
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz, :));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);
uvp_thresh  = 0.01 * max(phi_uvp_ref);
mask_ok     = phi_uvp_ref >= uvp_thresh;

% ---------------------------------------------------------------
% 3. Durkin trap data
% ---------------------------------------------------------------
trap = load_durkin_flux(trap_file);
mask_trap = trap.d_um >= 100 & trap.d_um < 2000;
d_trap    = trap.d_um(mask_trap);
flux_trap_agg = trap.flux_agg(:, mask_trap);
flux_trap_fp  = trap.flux_fp(:, mask_trap);

trap_depths  = [125, 330, 500];
id_trap      = zeros(1, 3);
for i = 1:3
    [~, id_trap(i)] = min(abs(trap.depths - trap_depths(i)));
end
z_trap = trap.depths(id_trap);

% total trap flux at each comparison depth
trap_total_agg = sum(flux_trap_agg(id_trap, :), 2, 'omitnan');  % 3 x 1
trap_total_fp  = sum(flux_trap_fp(id_trap, :),  2, 'omitnan');

% ---------------------------------------------------------------
% 4. Run model (flux BC)
% ---------------------------------------------------------------
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc, :) = Y(k_bc, :) + flux_src;
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('Converged at cycle %d\n', icyc); break;
    end
end

% final run: accumulate profiles on cast days
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);
Y_sum   = zeros(n_z, cfg.n_sections);
Yfp_sum = zeros(n_z, cfg.n_sections);
phi_mod_uvp = zeros(numel(k_plot), 1);   % for UVP comparison
n_cast = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc, :) = Y(k_bc, :) + flux_src;
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Y_sum   = Y_sum   + Y;
        Yfp_sum = Yfp_sum + Yfp;
        Ytot_day = Y + Yfp;
        phi_mod_uvp = phi_mod_uvp + sum(Ytot_day(k_plot, mask_uvp), 2);
        n_cast = n_cast + 1;
    end
end
Y_mean   = Y_sum   / max(n_cast, 1);
Yfp_mean = Yfp_sum / max(n_cast, 1);
phi_mod_uvp_avg = phi_mod_uvp / max(n_cast, 1);

% model BV flux profile [BV m-2 d-1] at every layer — for Martin b
% (BV flux is analogous to POC flux; canonical b=0.858 is from POC)
F_bv_full = zeros(n_z, 1);
for k = 1:n_z
    F_bv_full(k) = sum(w_bin .* Y_mean(k, :));
end

% model number flux [particles m-2 d-1], 100-2000 um only — for trap comparison
% F_n = w * Y / V,  V = pi/6 * d^3
F_n_mod  = zeros(n_z, 1);   % aggregate
F_n_fp   = zeros(n_z, 1);   % fecal
for k = 1:n_z
    F_n_mod(k) = sum(w_bin(mask_uvp) .* Y_mean(k, mask_uvp)   ./ V_bin(mask_uvp));
    F_n_fp(k)  = sum(w_bin(mask_uvp) .* Yfp_mean(k, mask_uvp) ./ V_bin(mask_uvp));
end

% UVP standing stock ratio
uvp_ratio = NaN(numel(k_plot), 1);
uvp_ratio(mask_ok) = phi_mod_uvp_avg(mask_ok) ./ phi_uvp_ref(mask_ok);

% model flux at trap comparison depths (number flux, 100-2000 um)
k_trap = zeros(1, 3);
for i = 1:3
    [~, k_trap(i)] = min(abs(z_centers - z_trap(i)));
end
mod_flux_at_trap = F_n_mod(k_trap);   % 3 x 1
fp_flux_at_trap  = F_n_fp(k_trap);

% Martin reference: use BV flux normalized to model at 125 m
[~, k125] = min(abs(z_centers - 125));
F_ref    = F_bv_full(k125);
b_martin = 0.858;
z_martin = z_centers(k125:end);
F_martin = F_ref * (z_martin / z_centers(k125)) .^ (-b_martin);

% model Martin b (100-500 m), from BV flux — consistent with run_flux_profile_test.m
[~, k100] = min(abs(z_centers - 100));
[~, k500] = min(abs(z_centers - 500));
F_top = F_bv_full(k100);  F_bot = F_bv_full(k500);
if F_top > 0 && F_bot > 0
    b_model = -log(F_bot/F_top) / log(z_centers(k500)/z_centers(k100));
else
    b_model = NaN;
end

% ---------------------------------------------------------------
% 5. Print summary table
% ---------------------------------------------------------------
fprintf('\n=== THREE-WAY COMPARISON SUMMARY ===\n');
fprintf('\n-- UVP standing stock ratio (model/UVP, 100-2000 um) --\n');
for ki = 1:numel(k_plot)
    r = uvp_ratio(ki);
    if isnan(r)
        fprintf('  %4.0f m: [UVP sparse]\n', z_mod(ki));
    else
        fprintf('  %4.0f m: %.2f\n', z_mod(ki), r);
    end
end

fprintf('\n-- Flux profile: model Martin b (100-500 m) --\n');
fprintf('  Model b = %.3f   (canonical = 0.858)\n', b_model);

fprintf('\n-- Durkin trap comparison (aggregate, 100-2000 um) --\n');
fprintf('  %-8s  %-12s  %-12s  %-8s\n', 'Depth', 'Mod agg', 'Trap agg', 'Ratio');
for i = 1:3
    fprintf('  %4.0f m    %9.2e    %9.2e    %6.2f\n', ...
        z_trap(i), mod_flux_at_trap(i), trap_total_agg(i), ...
        mod_flux_at_trap(i) / max(trap_total_agg(i), 1e-30));
end

% ---------------------------------------------------------------
% 6. Figure: 2-panel
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 16 14], 'Color', 'white');

% --- Panel 1: flux profile ---
subplot(1, 2, 1);
hold on;

% Martin reference
semilogy(F_martin, z_martin, 'r:', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('Martin b=%.2f', b_martin));

% model aggregate flux (number flux, 100-2000 um)
semilogy(F_n_mod(k_plot), z_mod, 'k-', 'LineWidth', 1.3, ...
    'DisplayName', 'Model agg');

% model fecal flux
semilogy(F_n_fp(k_plot), z_mod, 'k--', 'LineWidth', 0.9, ...
    'DisplayName', 'Model fp');

% Durkin trap aggregate
semilogy(trap_total_agg, z_trap, 'bs', 'MarkerSize', 6, 'LineWidth', 1.2, ...
    'DisplayName', 'Trap agg');

% Durkin trap fecal
semilogy(trap_total_fp, z_trap, 'b^', 'MarkerSize', 5, 'LineWidth', 1.0, ...
    'DisplayName', 'Trap fp');

set(gca, 'YDir', 'reverse', 'YLim', [60 550], 'FontSize', 7);
xlabel('Flux (particles m^{-2} d^{-1})');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Flux profile', 'FontWeight', 'normal');
hold off;

% --- Panel 2: ratio vs depth ---
subplot(1, 2, 2);
hold on;

% UVP standing stock ratio
r_uvp = uvp_ratio;
r_uvp(~mask_ok) = NaN;
plot(r_uvp, z_mod, 'g-o', 'MarkerSize', 4, 'LineWidth', 1.2, ...
    'DisplayName', 'Model/UVP (stock)');

% Durkin aggregate ratio
r_trap = mod_flux_at_trap ./ max(trap_total_agg, 1e-30);
plot(r_trap, z_trap, 'b-s', 'MarkerSize', 5, 'LineWidth', 1.2, ...
    'DisplayName', 'Model/Trap (flux)');

% ratio = 1 line
plot([1 1], [60 550], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

set(gca, 'YDir', 'reverse', 'YLim', [60 550], 'XScale', 'log', ...
    'XLim', [0.005 20], 'FontSize', 7);
xlabel('Model / Data ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Model/Data ratio', 'FontWeight', 'normal');
hold off;

saveas(gcf, fullfile(fig_dir, '3way_comparison.png'));
fprintf('\nSaved 3way_comparison.png\n');

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

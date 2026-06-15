% run_mass_fraction_diagnostic.m
%
% At steady state, compute what fraction of model BV at each depth is in:
%   small  (<100 um)    : too small for UVP
%   mid    (100-2000 um): UVP-visible
%   large  (>2000 um)   : too large for UVP
%
% If the deep residual is caused by mass piling in large bins (disagg off),
% we will see frac_large increase with depth while frac_mid drops.
%
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

bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 3:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% model bin diameters
grid_c = cfg.derive();
av_vol = grid_c.av_vol(:);            % cm^3 per particle, [n_sec x 1]
r_cm   = (0.75 / pi * av_vol).^(1/3);
d_um   = 2 * r_cm * 1e4;             % diameter [um]

mask_small = d_um < 100;
mask_mid   = d_um >= 100 & d_um < 2000;
mask_large = d_um >= 2000;

fprintf('Bins:  small=%d (<%.0f um)  mid=%d (%.0f-%.0f um)  large=%d (>%.0f um)\n', ...
    sum(mask_small), max(d_um(mask_small)), ...
    sum(mask_mid),   min(d_um(mask_mid)), max(d_um(mask_mid)), ...
    sum(mask_large), min(d_um(mask_large)));

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
% 3. Final run — mean over cast days, all depths
% ---------------------------------------------------------------
Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);
Y_sum = zeros(col_grid.n_z, cfg.n_sections);
n_cast = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    for i_step = 1:steps_per_day
        Y(k_bc, :) = phi_bc_daily(i_day, :);
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        Y(k_bc, :) = phi_bc_daily(i_day, :);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Y_sum  = Y_sum + Y + Yfp;
        n_cast = n_cast + 1;
    end
end
Ymean = Y_sum / max(n_cast, 1);   % [n_z x n_sec], averaged over cast days
fprintf('Model run complete. Cast days: %d\n', n_cast);

% ---------------------------------------------------------------
% 4. Compute BV in each size class at each depth
% ---------------------------------------------------------------
z_all    = col_grid.z_centers;

BV_small = sum(Ymean(:, mask_small), 2);
BV_mid   = sum(Ymean(:, mask_mid),   2);
BV_large = sum(Ymean(:, mask_large), 2);
BV_total = BV_small + BV_mid + BV_large;

safe_tot   = max(BV_total, 1e-30);
frac_small = 100 * BV_small ./ safe_tot;
frac_mid   = 100 * BV_mid   ./ safe_tot;
frac_large = 100 * BV_large ./ safe_tot;

% ---------------------------------------------------------------
% 5. Print table
% ---------------------------------------------------------------
fprintf('\nDepth    <100um    100-2000um   >2000um    BV_total\n');
for k = 1:col_grid.n_z
    fprintf('%5.0f m   %5.1f%%      %5.1f%%       %5.1f%%     %.2e\n', ...
        z_all(k), frac_small(k), frac_mid(k), frac_large(k), BV_total(k));
end

% ---------------------------------------------------------------
% 6. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 14 10], 'Color', 'white');

ax1 = subplot(1, 2, 1);
plot(frac_small, z_all, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
plot(frac_mid,   z_all, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2);
plot(frac_large, z_all, 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2);
set(ax1, 'YDir', 'reverse', 'YLim', [25 1000], 'XLim', [0 100]);
xlabel('% of total BV');
ylabel('Depth (m)');
legend('<100 \mum', '100-2000 \mum', '>2000 \mum', 'Location', 'southeast', 'FontSize', 7);
title('a) BV fractions', 'FontWeight', 'normal');

ax2 = subplot(1, 2, 2);
semilogy(BV_mid   * 1e6, z_all, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
semilogy(BV_large * 1e6, z_all, 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2);
semilogy(BV_small * 1e6, z_all, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2);
set(ax2, 'YDir', 'reverse', 'YLim', [25 1000]);
xlabel('BV (ppmV)');
legend('100-2000 \mum', '>2000 \mum', '<100 \mum', 'Location', 'southeast', 'FontSize', 7);
title('b) BV profiles', 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'mass_fraction_diagnostic.png'));
fprintf('\nSaved mass_fraction_diagnostic.png\n');

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

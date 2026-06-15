% run_100m_bin_range_check.m
%
% Find where the missing deep mass is.
%
% After disagg-off test: disaggregation is keeping mass IN the UVP window,
% not removing it. Coagulation pushes mass above 2 mm even with disagg on.
%
% Key question: at each depth, how much model mass is in:
%   (a) <100 um  (below UVP detection)
%   (b) 100-2000 um  (UVP window -- what we compare)
%   (c) >2000 um  (above UVP -- too big to see)
%
% If (c) is large at depth, then coagulation is moving mass above the UVP
% window faster than disaggregation can recycle it back. The fix is then
% either lower alpha or a stronger D_max constraint at depth.
%
% Uses baseline config from run_100m_start (all physics on).

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% ---------------------------------------------------------------
% 1. Setup (same baseline as run_100m_start)
% ---------------------------------------------------------------
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);

k_bc      = 2;
dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

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
cfg.surface_pp_mu    = 0.0;
cfg.r_to_rg          = 1.6;
cfg.zoo_c            = 0.025;
cfg.zoo_s            = 1.3e-5;
cfg.zoo_p            = 0.5;
cfg.zoo_ic           = 7;
cfg.mining_s         = 1.3e-5;
cfg.fp_alpha_cross   = 0.5;
cfg.validate();

grid_cfg   = cfg.derive();
r_cm       = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;
n_sec      = cfg.n_sections;

% bin range masks
mask_small = d_model_um < 100;
mask_uvp   = d_model_um >= 100 & d_model_um < 2000;
mask_large = d_model_um >= 2000;
fprintf('Bins <100 um: %d,  100-2000 um: %d,  >2000 um: %d\n', ...
    sum(mask_small), sum(mask_uvp), sum(mask_large));

% ---------------------------------------------------------------
% 2. BC at 100 m (same as run_100m_start)
% ---------------------------------------------------------------
uvpd = parse_uvp_daily(uvp_file);
[~, iz_100] = min(abs(uvpd.depth_m - 100));

uvp_ok    = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_ok  = uvpd.d_um(uvp_ok);
dw_uvp_ok = uvpd.dw(uvp_ok);
n_uvp_ok  = sum(uvp_ok);

d_model_edges = zeros(1, n_sec + 1);
d_model_edges(1)       = d_model_um(1)^2 / d_model_um(2);
d_model_edges(n_sec+1) = d_model_um(n_sec)^2 / d_model_um(n_sec-1);
for k = 2:n_sec
    d_model_edges(k) = sqrt(d_model_um(k-1) * d_model_um(k));
end

overlap_frac = zeros(n_sec, n_uvp_ok);
for j = 1:n_uvp_ok
    uvp_lo = d_uvp_ok(j) - dw_uvp_ok(j);
    uvp_hi = d_uvp_ok(j);
    for k = 1:n_sec
        lo = max(d_model_edges(k),   uvp_lo);
        hi = min(d_model_edges(k+1), uvp_hi);
        if hi > lo
            overlap_frac(k,j) = (hi - lo) / dw_uvp_ok(j);
        end
    end
end

n_days_uvp   = numel(uvpd.dates);
phi_100m     = zeros(n_days_uvp, n_sec);
for id = 1:n_days_uvp
    phi_row = squeeze(uvpd.phi(id, iz_100, uvp_ok));
    phi_row(isnan(phi_row)) = 0;
    for k = 1:n_sec
        phi_100m(id, k) = sum(overlap_frac(k,:) .* phi_row(:)');
    end
end

daily_surf   = get_daily_surface_phi(uvp_file, cfg, col_grid);
n_days       = daily_surf.n_days;
phi_bc_daily = zeros(n_days, n_sec);
for id_m = 1:n_days
    [~, ~, ib] = intersect(daily_surf.dates(id_m), uvpd.dates);
    if ~isempty(ib)
        phi_bc_daily(id_m, :) = phi_100m(ib, :);
    else
        dn_m   = datenum(num2str(daily_surf.dates(id_m)), 'yyyymmdd');
        dn_uvp = datenum(num2str(uvpd.dates), 'yyyymmdd');
        [~, nearest] = min(abs(dn_uvp - dn_m));
        phi_bc_daily(id_m, :) = phi_100m(nearest, :);
    end
end

[~, ia_match, ib_match] = intersect(daily_surf.dates, uvpd.dates);
phi_100m_total = sum(phi_bc_daily(ia_match, :), 2);
[~, ibest]     = max(phi_100m_total);
id_model_best  = ia_match(ibest);
id_uvp_best    = ib_match(ibest);
fprintf('Best cast day: %d\n', daily_surf.dates(id_model_best));

% ---------------------------------------------------------------
% 3. Run spinup + snapshot
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
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    if rel_change < spinup_tol
        fprintf('Converged at cycle %d\n', icyc);
        break;
    end
end

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
% 4. Bin-range breakdown at each depth
% ---------------------------------------------------------------
n_z   = col_grid.n_z;
k_cmp = 3:10;
z_cmp = col_grid.z_centers(k_cmp);

phi_small = sum(Y_snap(:, mask_small), 2);
phi_uvp_m = sum(Y_snap(:, mask_uvp),   2);
phi_large = sum(Y_snap(:, mask_large), 2);
phi_total = phi_small + phi_uvp_m + phi_large;

% UVP measured at comparison depths
phi_uvp_obs = zeros(numel(k_cmp), 1);
for i = 1:numel(k_cmp)
    [~, iz_u] = min(abs(uvpd.depth_m - z_cmp(i)));
    phi_row = squeeze(uvpd.phi(id_uvp_best, iz_u, uvp_ok));
    phi_row(isnan(phi_row)) = 0;
    phi_uvp_obs(i) = sum(phi_row);
end

fprintf('\n--- Bin-range breakdown by depth (best cast %d) ---\n', ...
    daily_surf.dates(id_model_best));
fprintf('%-8s  %-10s  %-10s  %-10s  %-10s  %-10s  %-8s\n', ...
    'z(m)', '<100um', '100-2000', '>2000um', 'total', 'UVP obs', 'ratio');
fprintf('%-8s  %-10s  %-10s  %-10s  %-10s  %-10s  %-8s\n', ...
    '', '[ppmV]', '[ppmV]', '[ppmV]', '[ppmV]', '[ppmV]', 'mod/UVP');
fprintf('%s\n', repmat('-', 1, 72));
for i = 1:numel(k_cmp)
    k = k_cmp(i);
    ratio = phi_uvp_m(k) / max(phi_uvp_obs(i), 1e-30);
    fprintf('%8.1f  %10.3f  %10.3f  %10.3f  %10.3f  %10.3f  %8.2f\n', ...
        z_cmp(i), ...
        phi_small(k)*1e6, phi_uvp_m(k)*1e6, phi_large(k)*1e6, ...
        phi_total(k)*1e6, phi_uvp_obs(i)*1e6, ratio);
end

% ---------------------------------------------------------------
% 5. Figure: three phi profiles stacked
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

figure('Units','centimeters','Position',[2 2 10 13]);
hold on;
plot(phi_small(k_cmp)*1e6,  z_cmp, 'b-',  'LineWidth', 1.5, 'DisplayName', '<100 \mum');
plot(phi_uvp_m(k_cmp)*1e6,  z_cmp, 'r-',  'LineWidth', 1.5, 'DisplayName', '100-2000 \mum (model)');
plot(phi_large(k_cmp)*1e6,  z_cmp, 'g-',  'LineWidth', 1.5, 'DisplayName', '>2000 \mum');
plot(phi_uvp_obs*1e6,        z_cmp, 'k--', 'LineWidth', 1.5, 'DisplayName', 'UVP obs');
set(gca, 'YDir', 'reverse');
xlabel('\phi [ppmV]');
ylabel('Depth (m)');
ylim([100 500]);
legend('location', 'southeast', 'FontSize', 7);
title('Bin-range breakdown: 100 m start');

saveas(gcf, fullfile(fig_dir, 'bin_range_breakdown.png'));
fprintf('\nSaved bin_range_breakdown.png\n');

% run_100m_loss_toggle.m
%
% Find which loss term drives the deep (>250 m) mass deficit.
%
% Setup (Adrian, June 11): start at 100 m, compare to UVP at 125-475 m.
% The 100m start run gave ratios 0.86-0.92 at 125-225 m but only 0.25 at
% 475 m. Something removes too much mass below 250 m.
%
% Strategy: run 4 cases.
%   case 1 (baseline): all physics on  (from run_100m_start)
%   case 2: zoo off   (enable_zoo = false)
%   case 3: microbe off (enable_microbe = false)
%   case 4: mining off  (enable_mining = false)
%
% Compare ratio model/UVP at each depth for each case.
% The case that flips deep ratios toward 1.0 is the dominant loss term.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% ---------------------------------------------------------------
% 1. Shared setup (same as run_100m_start)
% ---------------------------------------------------------------
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);

k_bc      = 2;        % BC at layer 2 (z = 75 m, closest to 100 m)
k_compare = 3:10;     % compare layers 125-475 m
z_compare = col_grid.z_centers(k_compare);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

% ---------------------------------------------------------------
% 2. Shared BC: UVP at 100 m mapped to model bins
% ---------------------------------------------------------------
uvpd = parse_uvp_daily(uvp_file);
[~, iz_100] = min(abs(uvpd.depth_m - 100));

% base config (all physics on)
cfg_base = SimulationConfig();
cfg_base.n_sections       = 30;
cfg_base.sinking_law      = 'kriest_8';
cfg_base.disagg_mode      = 'operator_split';
cfg_base.disagg_dmax_cm   = 1.0;
cfg_base.disagg_dmax_A    = 9.39e-6 * 5;
cfg_base.enable_coag      = true;
cfg_base.enable_disagg    = true;
cfg_base.enable_zoo       = true;
cfg_base.enable_microbe   = true;
cfg_base.enable_mining    = true;
cfg_base.alpha            = 0.5;
cfg_base.microbe_r0       = 0.03;
cfg_base.microbe_use_temp = true;
cfg_base.microbe_tref_C   = 20;
cfg_base.surface_pp_mu    = 0.0;
cfg_base.r_to_rg          = 1.6;
cfg_base.zoo_c            = 0.025;
cfg_base.zoo_s            = 1.3e-5;
cfg_base.zoo_p            = 0.5;
cfg_base.zoo_ic           = 7;
cfg_base.mining_s         = 1.3e-5;
cfg_base.fp_alpha_cross   = 0.5;
cfg_base.validate();

% model bin diameters (um)
grid_cfg   = cfg_base.derive();
r_cm       = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;
n_sec      = cfg_base.n_sections;

% overlap fraction (UVP 100-2000 um -> model bins)
uvp_ok   = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_ok = uvpd.d_um(uvp_ok);
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

% 100m phi in model bins for each UVP cast day
n_days_uvp = numel(uvpd.dates);
phi_100m = zeros(n_days_uvp, n_sec);
for id = 1:n_days_uvp
    phi_row = squeeze(uvpd.phi(id, iz_100, uvp_ok));
    phi_row(isnan(phi_row)) = 0;
    for k = 1:n_sec
        phi_100m(id, k) = sum(overlap_frac(k,:) .* phi_row(:)');
    end
end

% daily BC time series (match model days to UVP cast dates)
daily_surf = get_daily_surface_phi(uvp_file, cfg_base, col_grid);
n_days = daily_surf.n_days;
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

% best cast day
[~, ia_match, ib_match] = intersect(daily_surf.dates, uvpd.dates);
phi_100m_total = sum(phi_bc_daily(ia_match, :), 2);
[~, ibest] = max(phi_100m_total);
id_model_best = ia_match(ibest);
id_uvp_best   = ib_match(ibest);
fprintf('Best cast day: %d\n', daily_surf.dates(id_model_best));

% UVP phi at comparison depths
phi_uvp_cmp = zeros(numel(k_compare), 1);
for i = 1:numel(k_compare)
    [~, iz_u] = min(abs(uvpd.depth_m - z_compare(i)));
    phi_row = squeeze(uvpd.phi(id_uvp_best, iz_u, uvp_ok));
    phi_row(isnan(phi_row)) = 0;
    phi_uvp_cmp(i) = sum(phi_row);
end

mask_uvp_model = d_model_um >= 100 & d_model_um < 2000;

% ---------------------------------------------------------------
% 3. Four toggle cases
% ---------------------------------------------------------------
cases = struct();
cases(1).label      = 'all on (baseline)';
cases(1).zoo        = true;
cases(1).microbe    = true;
cases(1).mining     = true;

cases(2).label      = 'zoo OFF';
cases(2).zoo        = false;
cases(2).microbe    = true;
cases(2).mining     = true;

cases(3).label      = 'microbe OFF';
cases(3).zoo        = true;
cases(3).microbe    = false;
cases(3).mining     = true;

cases(4).label      = 'mining OFF';
cases(4).zoo        = true;
cases(4).microbe    = true;
cases(4).mining     = false;

n_cases = numel(cases);
ratio_table = NaN(numel(k_compare), n_cases);

for ic = 1:n_cases
    fprintf('\n=== Case %d: %s ===\n', ic, cases(ic).label);

    cfg = cfg_base;
    cfg.enable_zoo     = cases(ic).zoo;
    cfg.enable_microbe = cases(ic).microbe;
    cfg.enable_mining  = cases(ic).mining;
    cfg.validate();

    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(col_grid.n_z, n_sec);
    Yfp = zeros(col_grid.n_z, n_sec);

    % spinup
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
            fprintf('  Converged at cycle %d\n', icyc);
            break;
        end
    end

    % snapshot pass
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

    phi_mod = sum(Y_snap(k_compare, mask_uvp_model), 2);
    ratio_table(:, ic) = phi_mod ./ max(phi_uvp_cmp, 1e-30);
end

% ---------------------------------------------------------------
% 4. Print ratio table
% ---------------------------------------------------------------
fprintf('\n--- Ratio table: model (UVP range) / UVP ---\n');
fprintf('%-10s', 'Depth(m)');
for ic = 1:n_cases
    fprintf('  %-18s', cases(ic).label);
end
fprintf('\n%s\n', repmat('-', 1, 10 + n_cases*20));
for i = 1:numel(k_compare)
    fprintf('%-10.0f', z_compare(i));
    for ic = 1:n_cases
        fprintf('  %-18.2f', ratio_table(i, ic));
    end
    fprintf('\n');
end

% ---------------------------------------------------------------
% 5. Figure: ratio vs depth for all cases
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
figure('Units','centimeters','Position',[2 2 10 13]);
hold on;
for ic = 1:n_cases
    plot(ratio_table(:, ic), z_compare, '-o', ...
        'Color', colors(ic,:), 'LineWidth', 1.5, ...
        'DisplayName', cases(ic).label);
end
xline(1.0, 'k--', 'perfect match');
set(gca, 'YDir', 'reverse');
xlabel('model / UVP');
ylabel('Depth (m)');
ylim([100 500]);
xlim([0 4]);
legend('location', 'northeast', 'FontSize', 7);
title('Loss term toggle: 100 m start');

saveas(gcf, fullfile(fig_dir, 'loss_toggle_ratio.png'));
fprintf('\nSaved loss_toggle_ratio.png\n');

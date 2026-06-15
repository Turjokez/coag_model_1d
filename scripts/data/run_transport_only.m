% run_transport_only.m
%
% Test whether pure sinking (transport only) explains the deep profile.
%
% From report_june12 Section 8: the missing deep mass is not hiding in any
% size range -- it has genuinely left the column. Suspect: sinking speed
% too fast, particles exporting out of each layer before accumulating.
%
% Three cases:
%   1. baseline (all physics on)
%   2. transport only (coag/disagg/zoo/microbe/mining/pp all off)
%   3. transport + disagg only (isolate the recycling pump effect)
%
% If transport-only gives LESS than UVP: sinking is too fast on its own.
% If transport-only gives SIMILAR to UVP: physics terms are hurting, not helping.
% If transport-only gives MORE than UVP: physics removes mass from the UVP window.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% ---------------------------------------------------------------
% 1. Shared setup
% ---------------------------------------------------------------
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);

k_bc      = 2;
k_compare = 3:10;
z_compare = col_grid.z_centers(k_compare);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

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

grid_cfg   = cfg_base.derive();
r_cm       = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;
n_sec      = cfg_base.n_sections;

% ---------------------------------------------------------------
% 2. BC at 100 m
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

daily_surf   = get_daily_surface_phi(uvp_file, cfg_base, col_grid);
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

phi_uvp_cmp = zeros(numel(k_compare), 1);
for i = 1:numel(k_compare)
    [~, iz_u] = min(abs(uvpd.depth_m - z_compare(i)));
    phi_row = squeeze(uvpd.phi(id_uvp_best, iz_u, uvp_ok));
    phi_row(isnan(phi_row)) = 0;
    phi_uvp_cmp(i) = sum(phi_row);
end

mask_uvp_model = d_model_um >= 100 & d_model_um < 2000;

% ---------------------------------------------------------------
% 3. Three cases
% ---------------------------------------------------------------
cases = struct();
cases(1).label   = 'baseline (all on)';
cases(1).coag    = true;
cases(1).disagg  = true;
cases(1).zoo     = true;
cases(1).microbe = true;
cases(1).mining  = true;

cases(2).label   = 'transport only';
cases(2).coag    = false;
cases(2).disagg  = false;
cases(2).zoo     = false;
cases(2).microbe = false;
cases(2).mining  = false;

cases(3).label   = 'transport + disagg';
cases(3).coag    = false;
cases(3).disagg  = true;
cases(3).zoo     = false;
cases(3).microbe = false;
cases(3).mining  = false;

n_cases     = numel(cases);
ratio_table = NaN(numel(k_compare), n_cases);
Y_snaps     = cell(n_cases, 1);

for ic = 1:n_cases
    fprintf('\n=== Case %d: %s ===\n', ic, cases(ic).label);

    cfg = cfg_base;
    cfg.enable_coag    = cases(ic).coag;
    cfg.enable_disagg  = cases(ic).disagg;
    cfg.enable_zoo     = cases(ic).zoo;
    cfg.enable_microbe = cases(ic).microbe;
    cfg.enable_mining  = cases(ic).mining;
    cfg.validate();

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
            fprintf('  Converged at cycle %d\n', icyc);
            break;
        end
    end

    Y   = zeros(col_grid.n_z, n_sec);
    Yfp = zeros(col_grid.n_z, n_sec);
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
        if i_day == id_model_best
            Y_snaps{ic} = Y + Yfp;
        end
    end

    phi_mod = sum(Y_snaps{ic}(k_compare, mask_uvp_model), 2);
    ratio_table(:, ic) = phi_mod ./ max(phi_uvp_cmp, 1e-30);
end

% ---------------------------------------------------------------
% 4. Print ratio table
% ---------------------------------------------------------------
fprintf('\n--- Ratio table: model (UVP range) / UVP ---\n');
fprintf('%-10s', 'Depth(m)');
for ic = 1:n_cases
    fprintf('  %-24s', cases(ic).label);
end
fprintf('\n%s\n', repmat('-', 1, 10 + n_cases*26));
for i = 1:numel(k_compare)
    fprintf('%-10.0f', z_compare(i));
    for ic = 1:n_cases
        fprintf('  %-24.2f', ratio_table(i, ic));
    end
    fprintf('\n');
end

% ---------------------------------------------------------------
% 5. Figure: phi profile
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
figure('Units','centimeters','Position',[2 2 9 13]);
hold on;
plot(phi_uvp_cmp * 1e6, z_compare, 'k--', 'LineWidth', 1.5, 'DisplayName', 'UVP');
for ic = 1:n_cases
    phi_mod = sum(Y_snaps{ic}(k_compare, mask_uvp_model), 2);
    plot(phi_mod * 1e6, z_compare, '-o', 'Color', colors(ic,:), ...
        'LineWidth', 1.5, 'DisplayName', cases(ic).label);
end
set(gca, 'YDir', 'reverse');
xlabel('\phi [ppmV]');
ylabel('Depth (m)');
ylim([100 500]);
legend('location', 'southeast', 'FontSize', 7);
title('Transport only: 100 m start');

saveas(gcf, fullfile(fig_dir, 'transport_only_profile.png'));
fprintf('\nSaved transport_only_profile.png\n');

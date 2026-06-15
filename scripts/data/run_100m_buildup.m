% run_100m_buildup.m
%
% Build up full physics one term at a time on top of transport.
% Find which step causes the deep ratio to drop sharply.
%
% From report_june12 Section 9: transport alone gives ratio ~1 at 125-375 m.
% Full physics gives 0.25-0.92. The extra loss comes from coag + bio coupling,
% not from any single term alone. Here we add processes one by one.
%
% Cases (cumulative):
%   1. transport only
%   2. + coagulation
%   3. + coagulation + microbe
%   4. + coagulation + zoo
%   5. + coagulation + microbe + zoo  (= full without mining)
%   6. + all  (baseline, add mining)
%
% The step where the ratio drops most = the dominant coupling.

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
% 3. Build-up cases
% ---------------------------------------------------------------
cases = struct();

cases(1).label   = 'transport only';
cases(1).coag    = false;
cases(1).disagg  = false;
cases(1).microbe = false;
cases(1).zoo     = false;
cases(1).mining  = false;

cases(2).label   = '+ coag';
cases(2).coag    = true;
cases(2).disagg  = true;   % disagg always on with coag (recycling pump)
cases(2).microbe = false;
cases(2).zoo     = false;
cases(2).mining  = false;

cases(3).label   = '+ coag + microbe';
cases(3).coag    = true;
cases(3).disagg  = true;
cases(3).microbe = true;
cases(3).zoo     = false;
cases(3).mining  = false;

cases(4).label   = '+ coag + zoo';
cases(4).coag    = true;
cases(4).disagg  = true;
cases(4).microbe = false;
cases(4).zoo     = true;
cases(4).mining  = false;

cases(5).label   = '+ coag + microbe + zoo';
cases(5).coag    = true;
cases(5).disagg  = true;
cases(5).microbe = true;
cases(5).zoo     = true;
cases(5).mining  = false;

cases(6).label   = 'all on (baseline)';
cases(6).coag    = true;
cases(6).disagg  = true;
cases(6).microbe = true;
cases(6).zoo     = true;
cases(6).mining  = true;

n_cases     = numel(cases);
ratio_table = NaN(numel(k_compare), n_cases);
Y_snaps     = cell(n_cases, 1);

for ic = 1:n_cases
    fprintf('\n=== Case %d: %s ===\n', ic, cases(ic).label);

    cfg = cfg_base;
    cfg.enable_coag    = cases(ic).coag;
    cfg.enable_disagg  = cases(ic).disagg;
    cfg.enable_microbe = cases(ic).microbe;
    cfg.enable_zoo     = cases(ic).zoo;
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
fprintf('\n--- Build-up ratio table: model (UVP range) / UVP ---\n');
fprintf('%-10s', 'Depth(m)');
for ic = 1:n_cases
    lbl = cases(ic).label;
    if numel(lbl) > 22, lbl = lbl(1:22); end
    fprintf('  %-22s', lbl);
end
fprintf('\n%s\n', repmat('-', 1, 10 + n_cases*24));
for i = 1:numel(k_compare)
    fprintf('%-10.0f', z_compare(i));
    for ic = 1:n_cases
        fprintf('  %-22.2f', ratio_table(i, ic));
    end
    fprintf('\n');
end

% ---------------------------------------------------------------
% 5. Figure: ratio vs depth for all cases
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

colors = lines(n_cases);
figure('Units','centimeters','Position',[2 2 11 13]);
hold on;
for ic = 1:n_cases
    plot(ratio_table(:, ic), z_compare, '-o', ...
        'Color', colors(ic,:), 'LineWidth', 1.5, ...
        'DisplayName', cases(ic).label);
end
xline(1.0, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('model / UVP');
ylabel('Depth (m)');
ylim([100 500]);
xlim([0 2.5]);
legend('location', 'northeast', 'FontSize', 6);
title('Physics build-up: 100 m start');

saveas(gcf, fullfile(fig_dir, 'buildup_ratio.png'));
fprintf('\nSaved buildup_ratio.png\n');

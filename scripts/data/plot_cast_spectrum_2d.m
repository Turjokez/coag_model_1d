% plot_cast_spectrum_2d.m
%
% Cast-by-cast 2D volume spectrum comparison: UVP vs model.
% Recreates the panel figure style from Adrian's June meeting.
%
% Panel a: UVP particle volume spectrum [ppmV mm^-1], one column per cast date.
% Panel b: Model best config (alpha=0.10, Da x5, 100m BC).
%
% x: ESD (mm), log scale.  y: depth (m).  color: log10(spec).
% Only days with actual UVP casts are shown.
%
% Steps:
%   1. Load config + BC at 100 m.
%   2. Spinup model + save daily snapshots.
%   3. Build UVP and model spectral density images per cast day.
%   4. Plot 2 x n_cast panel figure with shared colorbar.

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
% 1. Config + grid
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
cfg.alpha          = 0.10;   % best from 2D grid
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
r_cm   = (0.75 / pi * grid_c.av_vol(:)).^(1/3);

% model bin diameters and widths [mm]
d_mod_mm = 2 * r_cm * 10;   % 1 cm = 10 mm
d_edges  = zeros(n_sec + 1, 1);
d_edges(1)       = d_mod_mm(1)^2   / d_mod_mm(2);
d_edges(n_sec+1) = d_mod_mm(n_sec)^2 / d_mod_mm(n_sec-1);
for k = 2:n_sec
    d_edges(k) = sqrt(d_mod_mm(k-1) * d_mod_mm(k));
end
dw_mod_mm = diff(d_edges);   % [n_sec x 1]

% ---------------------------------------------------------------
% 2. BC at 100 m + UVP struct
% ---------------------------------------------------------------
bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 3:10);
phi_bc_daily = bc.phi_bc_daily;   % [n_days x n_sec]
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% UVP bins 100-2000 um
mask_uvp  = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_mm  = uvpd.d_um(mask_uvp)  / 1000;   % [mm]
dw_uvp_mm = uvpd.dw(mask_uvp)    / 1000;   % [mm]

% UVP depth rows that fall inside our y plot range
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
% 4. Final run: save daily snapshots at model plot layers
% ---------------------------------------------------------------
Y   = zeros(col_grid.n_z, n_sec);
Yfp = zeros(col_grid.n_z, n_sec);
Y_daily = zeros(numel(k_plot), n_sec, n_days);   % [n_z x n_sec x n_days]

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
% 5. Find days with actual UVP casts
% ---------------------------------------------------------------
dn_obs = datenum(num2str(uvpd.dates), 'yyyymmdd');
dn_all = (dn_obs(1):dn_obs(end))';
[~, ia, ib] = intersect(bc.dates, uvpd.dates);
n_cast = numel(ia);
fprintf('Cast days with UVP data: %d\n', n_cast);

% date labels mm-dd
ds_all    = num2str(bc.dates(ia));
date_lbls = cell(n_cast, 1);
for m = 1:n_cast
    date_lbls{m} = [ds_all(m, 5:6) '-' ds_all(m, 7:8)];
end

% ---------------------------------------------------------------
% 6. Figure: 2 rows x n_cast columns
% ---------------------------------------------------------------
% color limits: log10(ppmV mm^-1)
cmin = -1;   cmax = 1;

% shared x limits (log10 mm)
xlim_log = [log10(0.09) log10(2.2)];
xtick_pos = log10([0.1 0.3 1.0]);
xtick_lbl = {'0.1', '0.3', '1.0'};

fig_w = max(n_cast * 1.3 + 1.5, 16);
figure('Units', 'centimeters', 'Position', [1 1 fig_w 10], 'Color', 'white');
colormap(jet);

% x coords in log10(mm) for each instrument
x_uvp_log = log10(d_uvp_mm);
x_mod_log  = log10(d_mod_mm);

for m = 1:n_cast
    id_mod = ia(m);   % index into bc.dates (model days 1..n_days)
    id_uvp = ib(m);   % index into uvpd.dates

    % --- UVP spectral density: [n_uvp_z x n_uvp_bins] ---
    phi_u = squeeze(uvpd.phi(id_uvp, mask_z_uvp, mask_uvp));   % [n_uvp_z x n_uvp_bins]
    phi_u(isnan(phi_u)) = 0;
    % parse_uvp_daily returns phi in [cm^3/cm^3], so:
    % S [ppmV mm^-1] = phi / dw_mm * 1e6
    if size(phi_u, 2) == numel(dw_uvp_mm)
        S_u = bsxfun(@rdivide, phi_u, dw_uvp_mm) * 1e6;
    elseif size(phi_u, 1) == numel(dw_uvp_mm)
        S_u = bsxfun(@rdivide, phi_u', dw_uvp_mm) * 1e6;
    else
        error('UVP spectrum shape does not match UVP bin widths.');
    end
    S_u(S_u <= 0) = NaN;

    % --- model spectral density: [n_z x n_sec] ---
    phi_m = squeeze(Y_daily(:, :, id_mod));   % [n_z x n_sec]
    S_m   = bsxfun(@rdivide, phi_m, dw_mod_mm') * 1e6;
    S_m(S_m <= 0) = NaN;

    % --- UVP panel (row 1) ---
    ax1 = subplot(2, n_cast, m);
    hU = imagesc(x_uvp_log, z_uvp, log10(S_u));
    set(hU, 'AlphaData', double(~isnan(S_u)));
    set(ax1, 'YDir', 'normal', 'CLim', [cmin cmax], ...
        'FontSize', 5, 'Color', 'white', ...
        'XLim', xlim_log, 'YLim', [60 510], ...
        'XTick', xtick_pos, 'XTickLabel', {});
    title(date_lbls{m}, 'FontSize', 5, 'FontWeight', 'normal');
    if m == 1
        ylabel('Depth (m)', 'FontSize', 7);
    else
        set(ax1, 'YTickLabel', {});
    end

    % --- Model panel (row 2) ---
    ax2 = subplot(2, n_cast, n_cast + m);
    hM = imagesc(x_mod_log, z_mod, log10(S_m));
    set(hM, 'AlphaData', double(~isnan(S_m)));
    set(ax2, 'YDir', 'normal', 'CLim', [cmin cmax], ...
        'FontSize', 5, 'Color', 'white', ...
        'XLim', xlim_log, 'YLim', [60 510], ...
        'XTick', xtick_pos, 'XTickLabel', xtick_lbl);
    if m == 1
        ylabel('Depth (m)', 'FontSize', 7);
        xlabel('ESD (mm)', 'FontSize', 6);
    else
        set(ax2, 'YTickLabel', {});
    end
end

% row labels
annotation('textbox', [0.005 0.82 0.03 0.08], 'String', 'a)', ...
    'EdgeColor', 'none', 'FontSize', 9, 'FontWeight', 'bold');
annotation('textbox', [0.005 0.35 0.03 0.08], 'String', 'b)', ...
    'EdgeColor', 'none', 'FontSize', 9, 'FontWeight', 'bold');

% colorbar
cb = colorbar('Position', [0.945 0.08 0.012 0.84]);
cb.Label.String = 'Particle Volume Spectrum (ppmV mm^{-1})';
cb.Label.FontSize = 6;
cb.Label.Interpreter = 'tex';
set(cb, 'Ticks', [-1 0 1], 'TickLabelInterpreter', 'tex', 'FontSize', 6);
cb.TickLabels = {'10^{-1}', '10^0', '10^1'};

saveas(gcf, fullfile(fig_dir, 'cast_spectrum_2d.png'));
fprintf('Saved cast_spectrum_2d.png\n');

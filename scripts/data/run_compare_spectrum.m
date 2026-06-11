% run_compare_spectrum.m
%
% Day-by-day, bin-by-bin size spectrum comparison: model vs UVP.
%
% For each cast day with UVP data, and at selected depths, plot the full
% particle size distribution (phi vs diameter) from model and UVP.
%
% Steps:
%   1. Run model with spinup (same config as run_data_column_daily)
%   2. Parse UVP by date and depth
%   3. Map model bins -> UVP bin space (integrate model phi within each UVP bin)
%   4. Figure A: spectrum at selected depth for all cast days (overlaid)
%   5. Figure B: ratio (model/UVP total phi) vs depth, per cast day
%   6. Figure C: one example day -- spectra at three depths

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- model config (same as run_data_column_daily) ---
cfg = SimulationConfig();
cfg.n_sections    = 30;
cfg.sinking_law   = 'kriest_8';
cfg.disagg_mode   = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.enable_coag   = true;
cfg.enable_disagg = true;
cfg.enable_zoo    = true;
cfg.enable_microbe   = true;
cfg.enable_mining    = true;
cfg.alpha            = 0.5;
cfg.microbe_r0       = 0.03;
cfg.microbe_use_temp = true;    % Q10=2 scaling with real T(z)
cfg.microbe_tref_C   = 20;      % reference temp [C]
cfg.surface_pp_mu  = 0.1;
cfg.r_to_rg        = 1.6;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.mining_s       = 1.3e-5;
cfg.fp_alpha_cross = 0.5;
cfg.disagg_dmax_A  = 9.39e-6 * 5;   % x5: D_max ~ 1.48 mm at surface eps
cfg.validate();

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

col_grid  = ColumnGrid(1000, 20);
prof      = load_keps(mat_path, col_grid.z_centers);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);   % daily eps
n_z       = col_grid.n_z;
n_sec     = cfg.n_sections;

% model bin diameters [um]
grid_cfg   = cfg.derive();
r_cm       = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_um = (2 * r_cm * 1e4)';   % 1 x n_sec

% model bin boundaries (geometric mean of adjacent centers, log-spaced)
log_d   = log(d_model_um);
log_bnd = [log_d(1) - (log_d(2)-log_d(1))/2, ...
           (log_d(1:end-1) + log_d(2:end)) / 2, ...
           log_d(end) + (log_d(end)-log_d(end-1))/2];
d_model_bounds_um = exp(log_bnd);   % 1 x (n_sec+1)

% --- spinup + comparison run ---
daily = get_daily_surface_phi(uvp_file, cfg, col_grid);
n_days = daily.n_days;

sim = ColumnSimulation(cfg, col_grid, prof);

fprintf('Spinup...\n');
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);
for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 3), 2);
    for i_day = 1:n_days
        for i_step = 1:steps_per_day
            Y(1, :) = daily.phi(i_day, :);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(1, :) = daily.phi(i_day, :);
        end
    end
    phi_after = mean(sum(Y + Yfp, 3), 2);
    rel_change = max(abs(phi_after - phi_before) ./ max(phi_before, 1e-20));
    if rel_change < spinup_tol
        fprintf('  converged at cycle %d\n', icyc);
        break;
    end
end

fprintf('Comparison run...\n');
Y_daily   = zeros(n_days, n_z, n_sec);
Yfp_daily = zeros(n_days, n_z, n_sec);
for i_day = 1:n_days
    % update eps profile for this day if data is available
    d_today = daily.dates(i_day);
    ikd = find(keps_day.dates == d_today, 1);
    if ~isempty(ikd)
        sim.rhs.profile.eps = keps_day.eps(:, ikd);
    end
    for i_step = 1:steps_per_day
        Y(1, :) = daily.phi(i_day, :);
        [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
        Y(1, :) = daily.phi(i_day, :);
    end
    Y_daily(i_day, :, :)   = Y;
    Yfp_daily(i_day, :, :) = Yfp;
end
fprintf('Done.\n');

% total model phi per day/depth/bin (agg + fecal)
model_phi = Y_daily + Yfp_daily;   % n_days x n_z x n_sec

% --- parse UVP by date and depth ---
uvpd = parse_uvp_daily(uvp_file);

% filter UVP to < 2000 um (same as main script)
uvp_bin_mask = uvpd.d_um >= 100 & uvpd.d_um < 2000;
n_uvp_bins   = sum(uvp_bin_mask);
d_uvp        = uvpd.d_um(uvp_bin_mask);

% --- map model bins -> UVP bins (weighted overlap) ---
% For each model bin k and UVP bin j, compute what fraction of model bin k
% falls inside UVP bin j. Then distribute phi proportionally.
% This avoids stripes from the old center-falls-in approach.
uvp_bounds = uvpd.d_bounds;   % 1 x (n_uvp_all+1) [um]
n_uvp_all  = numel(uvpd.d_um);

% overlap_frac(k, j) = fraction of model bin k that overlaps with UVP bin j
overlap_frac = zeros(n_sec, n_uvp_all);
for k = 1:n_sec
    lo_k = d_model_bounds_um(k);
    hi_k = d_model_bounds_um(k+1);
    for j = 1:n_uvp_all
        ov = max(0, min(hi_k, uvp_bounds(j+1)) - max(lo_k, uvp_bounds(j)));
        overlap_frac(k, j) = ov / (hi_k - lo_k);
    end
end

% project model phi into UVP bin space
% model_in_uvp(day, depth, uvp_bin) = sum_k phi_k * overlap_frac(k, uvp_bin)
model_in_uvp = zeros(n_days, n_z, n_uvp_all);
for k = 1:n_sec
    for j = 1:n_uvp_all
        if overlap_frac(k, j) > 0
            model_in_uvp(:, :, j) = model_in_uvp(:, :, j) + ...
                model_phi(:, :, k) * overlap_frac(k, j);
        end
    end
end

% --- match model days to UVP cast dates ---
% daily.dates is YYYYMMDD for model day indices 1..n_days
% uvpd.dates  is YYYYMMDD for cast days
[cast_dates_matched, ia, ib] = intersect(daily.dates, uvpd.dates);
n_matched = numel(cast_dates_matched);
fprintf('Matched %d model days to UVP casts.\n', n_matched);

% --- Figure A: spectrum at one depth for all matched cast days ---
% pick depth closest to 150 m
[~, iz_sel] = min(abs(col_grid.z_centers - 150));
z_sel = col_grid.z_centers(iz_sel);

figure;
hold on;
cmap = lines(n_matched);
for m = 1:n_matched
    id_model = ia(m);   % model day index
    id_uvp   = ib(m);   % uvpd.dates index

    % find nearest UVP depth to z_sel
    [~, iz_uvp] = min(abs(uvpd.depth_m - z_sel));

    phi_uvp   = squeeze(uvpd.phi(id_uvp, iz_uvp, uvp_bin_mask));
    phi_uvp(isnan(phi_uvp)) = 0;
    phi_model = squeeze(model_in_uvp(id_model, iz_sel, uvp_bin_mask));

    % only plot non-zero bins to keep log scale clean
    ok_m = phi_model > 0;
    ok_u = phi_uvp   > 0;
    if any(ok_m)
        loglog(d_uvp(ok_m), phi_model(ok_m), '-',  'Color', cmap(m, :), 'LineWidth', 1);
    end
    if any(ok_u)
        loglog(d_uvp(ok_u), phi_uvp(ok_u),   '--', 'Color', cmap(m, :), 'LineWidth', 1);
    end
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('diameter  [\mum]');
ylabel('\phi  [cm^3 cm^{-3}]');
title(sprintf('spectrum at z = %d m  (solid=model, dashed=UVP)', round(z_sel)));
saveas(gcf, fullfile(fig_dir, 'spectrum_at_150m.png'));

% --- Figure B: total phi ratio (model/UVP) vs depth, per cast day ---
% interpolate UVP phi to model z grid
figure;
hold on;
for m = 1:n_matched
    id_model = ia(m);
    id_uvp   = ib(m);

    % UVP total phi at each depth (filtered bins)
    uvp_phi_cast  = squeeze(uvpd.phi(id_uvp, :, uvp_bin_mask));  % n_uvp_depths x n_filtered_bins
    uvp_phi_cast(isnan(uvp_phi_cast)) = 0;
    uvp_phi_total = sum(uvp_phi_cast, 2);   % n_uvp_depths x 1

    % skip cast if fewer than 2 depths have any data
    valid = uvp_phi_total > 0;
    if sum(valid) < 2, continue; end

    % interpolate only within UVP cast depth range (no extrapolation)
    uvp_phi_interp = interp1(uvpd.depth_m(valid), uvp_phi_total(valid), ...
        col_grid.z_centers, 'pchip');   % NaN outside range
    uvp_phi_interp = max(0, uvp_phi_interp(:));   % column vector

    model_phi_total = squeeze(sum(model_in_uvp(id_model, :, uvp_bin_mask), 3));
    model_phi_total = model_phi_total(:);   % column vector

    % only compute ratio at depths with valid UVP data
    valid_z = ~isnan(uvp_phi_interp) & uvp_phi_interp > 1e-20;
    ratio = nan(n_z, 1);
    ratio(valid_z) = model_phi_total(valid_z) ./ uvp_phi_interp(valid_z);
    plot(ratio, col_grid.z_centers, '-', 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
end
xline(1, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('model / UVP  \phi ratio');
ylabel('depth  [m]');
title('model/UVP ratio per cast day (grey lines)');
saveas(gcf, fullfile(fig_dir, 'spectrum_ratio_vs_depth.png'));

% --- Figure C: one example cast day -- spectrum at three depths ---
% pick the cast day with highest surface forcing (most particles)
[~, best_cast] = max(sum(daily.phi(ia, :), 2));
id_model = ia(best_cast);
id_uvp   = ib(best_cast);
cast_date_str = num2str(cast_dates_matched(best_cast));

check_depths = [75, 200, 400];
styles = {'b-', 'r-', 'g-'};

figure;
for kd = 1:3
    [~, iz_uvp] = min(abs(uvpd.depth_m - check_depths(kd)));
    [~, iz_mod] = min(abs(col_grid.z_centers - check_depths(kd)));

    phi_uvp   = squeeze(uvpd.phi(id_uvp, iz_uvp, uvp_bin_mask));
    phi_uvp(isnan(phi_uvp)) = 0;
    phi_model = squeeze(model_in_uvp(id_model, iz_mod, uvp_bin_mask));

    loglog(d_uvp, phi_model, styles{kd}, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('model %dm', check_depths(kd)));
    hold on;
    loglog(d_uvp, phi_uvp, [styles{kd}(1) '--'], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('UVP %dm', check_depths(kd)));
end
xlabel('diameter  [\mum]');
ylabel('\phi  [cm^3 cm^{-3}]');
legend('location', 'southwest');
title(sprintf('spectra: model (solid) vs UVP (dashed)  date=%s', cast_date_str));
saveas(gcf, fullfile(fig_dir, 'spectrum_example_day.png'));

% --- Figure D: cast-by-cast spectrum panels (UVP top, model bottom) ---
% Units: ppmV/mm (differential particle volume spectrum, Siegel et al. 2025).
% Conversion: DVSD [uL/m3/um] = phi [cm3/cm3] / (dw [um] * 1e-9).
% Since 1 uL/m3/um = 1 ppmV/mm, S [ppmV/mm] = phi / (dw * 1e-9) = phi/dw * 1e9.
% x-axis: ESD in mm (same as reference figure).
% Fixed color range: 10^-1 to 10^1 ppmV/mm (same as Siegel et al.).
depth_lim = 500;

iz_mod = col_grid.z_centers <= depth_lim;
iz_uvp = uvpd.depth_m      <= depth_lim;
z_mod  = col_grid.z_centers(iz_mod);
z_uvp  = uvpd.depth_m(iz_uvp);

% bin widths for filtered bins [um], for spectrum conversion
dw_filt = uvpd.dw(uvp_bin_mask);   % 1 x n_uvp_bins

% x-axis in mm (log10 scale), matching Siegel et al.
d_uvp_mm = d_uvp / 1000;
x_log_mm = log10(d_uvp_mm);

% fixed color limits matching Siegel et al. Fig 2 colorbar: 10^-1 to 10^1 ppmV/mm
cmin = -1;
cmax =  1;

figure('Units', 'centimeters', 'Position', [2 2 max(16, n_matched*1.1) 8], ...
    'Color', 'white');
colormap(jet);

for m = 1:n_matched
    id_model = ia(m);
    id_uvp   = ib(m);
    ds  = num2str(cast_dates_matched(m));
    lbl = [ds(5:6) '-' ds(7:8)];

    % UVP: convert phi -> spectrum [ppmV/mm]
    phi_u = squeeze(uvpd.phi(id_uvp, iz_uvp, uvp_bin_mask));   % n_depths x n_bins
    S_u   = phi_u ./ dw_filt * 1e9;   % ppmV/mm, broadcast over bins
    S_u(S_u <= 0) = NaN;

    % model: convert phi -> spectrum [ppmV/mm]
    phi_m = squeeze(model_in_uvp(id_model, iz_mod, uvp_bin_mask));   % n_z x n_bins
    S_m   = phi_m ./ dw_filt * 1e9;
    S_m(S_m <= 0) = NaN;

    % top row: UVP
    ax = subplot(2, n_matched, m);
    h = imagesc(x_log_mm, z_uvp, log10(S_u));
    set(h, 'AlphaData', double(~isnan(S_u)));   % NaN -> transparent -> white bg
    set(ax, 'YDir', 'reverse', 'CLim', [cmin cmax], 'FontSize', 5, 'Color', 'white');
    title(lbl, 'FontSize', 6, 'FontWeight', 'normal');
    set(ax, 'XTickLabel', {});
    if m == 1
        ylabel('Depth (m)', 'FontSize', 7);
    else
        set(ax, 'YTickLabel', {});
    end

    % bottom row: model
    ax2 = subplot(2, n_matched, n_matched + m);
    h2 = imagesc(x_log_mm, z_mod, log10(S_m));
    set(h2, 'AlphaData', double(~isnan(S_m)));
    set(ax2, 'YDir', 'reverse', 'CLim', [cmin cmax], 'FontSize', 5, 'Color', 'white');
    set(ax2, 'XTick', log10([0.1 1]));
    if m == 1
        set(ax2, 'XTickLabel', {'10^{-1}', '10^0'}, 'TickLabelInterpreter', 'tex');
        xlabel('ESD (mm)', 'FontSize', 7);
        ylabel('Depth (m)', 'FontSize', 7);
    else
        set(ax2, 'XTickLabel', {}, 'YTickLabel', {});
    end
end

annotation('textbox', [0.01 0.88 0.04 0.08], 'String', 'a)', ...
    'EdgeColor', 'none', 'FontSize', 8, 'FontWeight', 'bold');
annotation('textbox', [0.01 0.40 0.04 0.08], 'String', 'b)', ...
    'EdgeColor', 'none', 'FontSize', 8, 'FontWeight', 'bold');

cb = colorbar('Position', [0.945 0.1 0.012 0.8]);
cb.Label.String = 'Particle Volume Spectrum (ppmV mm^{-1})';
cb.Label.FontSize = 6;
cb.Label.Interpreter = 'tex';
set(cb, 'Ticks', [-1 0 1], 'TickLabelInterpreter', 'tex');
cb.TickLabels = {'10^{-1}', '10^0', '10^1'};

saveas(gcf, fullfile(fig_dir, 'spectrum_cast_panels.png'));

fprintf('Figures saved to %s\n', fig_dir);

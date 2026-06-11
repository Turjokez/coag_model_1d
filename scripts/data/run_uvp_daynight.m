% run_uvp_daynight.m
%
% 3-row comparison: Night UVP / Day UVP / Model
% for 4 selected dates (May 7, 13, 19, 22).
% Day = UTC 06-21h, Night = UTC 21-06h (EXPORTS-NA, ~49N 16W).
% Bins filtered to 100-2000 um (UVP reliable range).

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', 'Turbulance', 'keps_for_dave.mat');
fig_dir  = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Read UVP file
% ---------------------------------------------------------------
fid = fopen(uvp_file, 'r');
missing_val = -9999;
fields_line = '';
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    if startsWith(line, '/missing='), missing_val = str2double(extractAfter(line,'/missing=')); end
    if startsWith(line, '/fields='),  fields_line = extractAfter(line, '/fields='); end
    if strcmp(strtrim(line), '/end_header'), break; end
end
all_fields = strtrim(strsplit(fields_line, ','));
col_depth = find(strcmp(all_fields, 'depth'));
col_date  = find(strcmp(all_fields, 'date'));
col_time  = find(strcmp(all_fields, 'time'));
dvsd_mask = strncmp(all_fields, 'PSD_DVSD_', 9);
col_dvsd  = find(dvsd_mask);
n_uvp     = numel(col_dvsd);

fmt = repmat('%s ', 1, numel(all_fields));
C   = textscan(fid, fmt, 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);
raw = C{1};

depth_all = cellfun(@str2double, raw(:, col_depth));
date_all  = cellfun(@str2double, raw(:, col_date));
hour_all  = cellfun(@(t) str2double(t(1:2)), raw(:, col_time));

DVSD = zeros(size(raw,1), n_uvp);
for j = 1:n_uvp
    vals = cellfun(@str2double, raw(:, col_dvsd(j)));
    vals(abs(vals - missing_val) < 1) = NaN;
    DVSD(:, j) = vals;
end

% ---------------------------------------------------------------
% 2. Bin geometry and UVP filter (100-2000 um)
% ---------------------------------------------------------------
d_bounds  = [50.8, 64, 80.6, 102, 128, 161, 203, 256, 323, 406, 512, 645, ...
             813, 1020, 1290, 1630, 2050, 2580, 3250, 4100, 5160, 6500, ...
             8190, 10300, 13000, 16400, 20600, 26000];
dw_um     = diff(d_bounds);                                 % bin widths, um
d_centers = sqrt(d_bounds(1:end-1) .* d_bounds(2:end));    % geometric centers, um
mask_uvp  = d_centers >= 100 & d_centers < 2000;
d_mm      = d_centers(mask_uvp) / 1000;                    % mm, for plotting
dw_filt   = dw_um(mask_uvp);                               % um

% depth grid for UVP averaging
z_edges = 0:10:500;
z_mid   = (z_edges(1:end-1) + z_edges(2:end)) / 2;
n_zdep  = numel(z_mid);
n_bins  = sum(mask_uvp);

% ---------------------------------------------------------------
% 3. Extract UVP day/night profiles for selected casts
% ---------------------------------------------------------------
sel      = [20210507, 1, 10;
            20210513, 2, 11;
            20210519, 0, 12;
            20210522, 4, 12];
n_dates  = size(sel, 1);
date_lbl = {'05-07', '05-13', '05-19', '05-22'};

S_uvp = NaN(n_dates, 2, n_zdep, n_bins);   % dim 2: 1=night, 2=day

for id = 1:n_dates
    d = sel(id, 1);
    for dn = 1:2
        target_h = sel(id, 1 + dn);   % night hour (col 2) or day hour (col 3)
        for iz = 1:n_zdep
            rows = date_all == d & hour_all == target_h & ...
                   depth_all >= z_edges(iz) & depth_all < z_edges(iz+1);
            if ~any(rows), continue; end
            dvsd_here = DVSD(rows, mask_uvp);
            dvsd_here(dvsd_here < 0) = NaN;
            S_uvp(id, dn, iz, :) = mean_no_nan(dvsd_here, 1);
        end
    end
end

% ---------------------------------------------------------------
% 4. Run model simulation (operator_split, daily eps)
% ---------------------------------------------------------------
cfg = SimulationConfig();
cfg.n_sections      = 30;
cfg.sinking_law     = 'kriest_8';
cfg.enable_coag     = true;
cfg.enable_disagg   = true;
cfg.disagg_mode     = 'operator_split';
cfg.disagg_dmax_cm  = 1.0;      % fallback; daily eps overwrites by depth/day
cfg.enable_zoo      = true;
cfg.enable_microbe  = true;
cfg.enable_mining   = true;
cfg.alpha           = 0.5;
cfg.microbe_r0      = 0.03;
cfg.microbe_use_temp = true;
cfg.microbe_tref_C  = 20;
cfg.surface_pp_mu   = 0.1;
cfg.r_to_rg         = 1.6;
cfg.zoo_c           = 0.025;
cfg.zoo_s           = 1.3e-5;
cfg.zoo_p           = 0.5;
cfg.zoo_ic          = 7;
cfg.mining_s        = 1.3e-5;
cfg.fp_alpha_cross  = 0.5;
cfg.validate();

col_grid  = ColumnGrid(1000, 20);
prof      = load_keps(mat_path, col_grid.z_centers);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);   % struct with .eps, .dates
daily     = get_daily_surface_phi(uvp_file, cfg, col_grid);
n_days    = daily.n_days;
n_z_mod   = col_grid.n_z;

% model -> UVP bin overlap fractions
grid_cfg = cfg.derive();
r_cm     = (0.75/pi * grid_cfg.av_vol(:)).^(1/3);
d_lo_um  = zeros(cfg.n_sections, 1);
d_hi_um  = zeros(cfg.n_sections, 1);
log_d    = log(2 * r_cm * 1e4);   % log(diameter in um)
log_bnd  = [log_d(1)-(log_d(2)-log_d(1))/2; ...
            (log_d(1:end-1)+log_d(2:end))/2; ...
            log_d(end)+(log_d(end)-log_d(end-1))/2];
d_lo_um  = exp(log_bnd(1:end-1));
d_hi_um  = exp(log_bnd(2:end));

n_uvp_all   = numel(d_centers);
overlap_frac = zeros(cfg.n_sections, n_uvp_all);
for k = 1:cfg.n_sections
    lo_k = d_lo_um(k); hi_k = d_hi_um(k);
    for j = 1:n_uvp_all
        lo_j = d_bounds(j); hi_j = d_bounds(j+1);
        ov = max(0, min(hi_k,hi_j) - max(lo_k,lo_j));
        overlap_frac(k, j) = ov / (hi_k - lo_k);
    end
end

% spinup
dt            = 0.25;
steps_per_day = round(1/dt);
spinup_tol    = 0.01;
max_cycles    = 50;

sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z_mod, cfg.n_sections);
Yfp = zeros(n_z_mod, cfg.n_sections);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(1,:) = daily.phi(i_day,:);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(1,:) = daily.phi(i_day,:);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    if rel_change < spinup_tol
        fprintf('Spinup converged at cycle %d\n', icyc);
        break;
    end
end

% extract snapshot for each selected date
sel_dates  = sel(:, 1);
S_model    = NaN(n_dates, n_zdep, n_bins);

% re-run one full pass to grab daily snapshots
Y   = zeros(n_z_mod, cfg.n_sections);
Yfp = zeros(n_z_mod, cfg.n_sections);
for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    for i_step = 1:steps_per_day
        Y(1,:) = daily.phi(i_day,:);
        [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
        Y(1,:) = daily.phi(i_day,:);
    end

    phi_snap = Y + Yfp;   % n_z_mod x n_sec

    for id = 1:n_dates
        if daily.dates(i_day) ~= sel_dates(id), continue; end

        % map to UVP bins: phi_uvp(n_z_mod, n_uvp_all)
        phi_uvp = zeros(n_z_mod, n_uvp_all);
        for k = 1:cfg.n_sections
            for j = 1:n_uvp_all
                if overlap_frac(k,j) > 0
                    phi_uvp(:, j) = phi_uvp(:, j) + phi_snap(:,k) * overlap_frac(k,j);
                end
            end
        end

        % interpolate model layers -> UVP depth grid (nearest neighbor)
        for iz = 1:n_zdep
            [~, kz] = min(abs(col_grid.z_centers - z_mid(iz)));
            phi_row = phi_uvp(kz, mask_uvp);           % 1 x n_bins
            S_model(id, iz, :) = phi_row ./ dw_filt * 1e9;   % ppmV/mm
        end
    end
end

% ---------------------------------------------------------------
% 5. Model with disagg OFF (4th row)
% ---------------------------------------------------------------
cfg_off = cfg.copy();
cfg_off.enable_disagg = false;
cfg_off.validate();

sim_off = ColumnSimulation(cfg_off, col_grid, prof);
Y   = zeros(n_z_mod, cfg.n_sections);
Yfp = zeros(n_z_mod, cfg.n_sections);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim_off.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(1,:) = daily.phi(i_day,:);
            [Y, Yfp] = sim_off.rhs.stepY(Y, dt, Yfp);
            Y(1,:) = daily.phi(i_day,:);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    if rel_change < spinup_tol
        fprintf('Disagg-off spinup converged at cycle %d\n', icyc);
        break;
    end
end

S_model_off = NaN(n_dates, n_zdep, n_bins);
Y   = zeros(n_z_mod, cfg.n_sections);
Yfp = zeros(n_z_mod, cfg.n_sections);
for i_day = 1:n_days
    sim_off.rhs.profile.eps = keps_day.eps(:, i_day);
    for i_step = 1:steps_per_day
        Y(1,:) = daily.phi(i_day,:);
        [Y, Yfp] = sim_off.rhs.stepY(Y, dt, Yfp);
        Y(1,:) = daily.phi(i_day,:);
    end
    phi_snap = Y + Yfp;
    for id = 1:n_dates
        if daily.dates(i_day) ~= sel_dates(id), continue; end
        phi_uvp = zeros(n_z_mod, n_uvp_all);
        for k = 1:cfg.n_sections
            for j = 1:n_uvp_all
                if overlap_frac(k,j) > 0
                    phi_uvp(:,j) = phi_uvp(:,j) + phi_snap(:,k) * overlap_frac(k,j);
                end
            end
        end
        for iz = 1:n_zdep
            [~, kz] = min(abs(col_grid.z_centers - z_mid(iz)));
            phi_row = phi_uvp(kz, mask_uvp);
            S_model_off(id, iz, :) = phi_row ./ dw_filt * 1e9;
        end
    end
end

% ---------------------------------------------------------------
% 6. Plot: 4 rows x 4 columns (Night / Day / Model / Model no disagg)
% ---------------------------------------------------------------
clim_log  = [-1 1];   % log10 ppmV/mm: 0.1 to 10
row_lbl   = {'Night', 'Day', 'Model', 'No disagg'};
n_rows    = 4;

figure('Units', 'centimeters', 'Position', [2 2 22 17]);

for row = 1:n_rows
    for id = 1:n_dates
        ax = subplot(n_rows, n_dates, (row-1)*n_dates + id);

        if row <= 2
            S = squeeze(S_uvp(id, row, :, :));
        elseif row == 3
            S = squeeze(S_model(id, :, :));
        else
            S = squeeze(S_model_off(id, :, :));
        end
        S(S <= 0) = NaN;
        Sl = log10(S);

        imagesc(log10(d_mm), z_mid, Sl);
        set(ax, 'YDir', 'reverse', 'XLim', [log10(0.095) 0.3]);
        clim(clim_log);
        colormap(ax, jet);

        x_ticks = log10([0.1 0.3 1]);
        set(ax, 'XTick', x_ticks);

        if row == n_rows
            set(ax, 'XTickLabel', {'0.1','0.3','1'});
            xlabel('ESD (mm)');
        else
            set(ax, 'XTickLabel', {});
        end
        if row == 1, title(date_lbl{id}); end

        if id == 1
            ylabel('Depth (m)');
            text(-0.3, 0.5, row_lbl{row}, 'Units','normalized', ...
                'Rotation', 90, 'HorizontalAlignment','center', 'FontSize', 8);
        else
            set(ax, 'YTickLabel', {});
        end
    end
end

cb = colorbar('Position', [0.92 0.1 0.02 0.8]);
cb.Label.String = 'S (ppmV mm^{-1})';
set(cb, 'Ticks', [-1 0 1], 'TickLabels', {'0.1','1','10'});

saveas(gcf, fullfile(fig_dir, 'uvp_daynight_model.png'));
fprintf('Saved uvp_daynight_model.png\n');

% ---------------------------------------------------------------
function y = mean_no_nan(x, dim)
if nargin < 2, dim = 1; end
mask = ~isnan(x);
num  = sum(x .* mask, dim);
den  = sum(mask, dim);
y = num ./ den;
y(den == 0) = NaN;
end

function daily = get_daily_surface_phi(sb_file, cfg, col_grid)
% GET_DAILY_SURFACE_PHI  Extract surface phi per day from UVP .sb file.
%
% Usage:
%   daily = get_daily_surface_phi(sb_file, cfg, col_grid)
%
% Output:
%   daily.day_num   - day index 1..n_days (from first date in file)
%   daily.dates     - YYYYMMDD for each day [n_days x 1]
%   daily.phi       - surface phi in model bins [n_days x n_sec]
%                     gaps filled by linear interpolation
%   daily.has_data  - logical [n_days x 1], true if real UVP data exists
%   daily.d_model   - model bin diameters [um]
%
% Notes:
%   - Only depth <= 5 m rows used (UVP 2.5 m bin).
%   - Multiple casts on same day are averaged.
%   - Missing days (no cast) filled by linear interpolation between
%     nearest neighbors. Edge gaps use nearest value.
%   - phi units: cm^3/cm^3 (same as model Y).

% --- read UVP file (surface rows only) ---
fid = fopen(sb_file, 'r');
if fid < 0, error('Cannot open: %s', sb_file); end

missing_val = -9999;
fields_line = '';
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    if startsWith(line, '/missing=')
        missing_val = str2double(extractAfter(line, '/missing='));
    end
    if startsWith(line, '/fields=')
        fields_line = extractAfter(line, '/fields=');
    end
    if strcmp(strtrim(line), '/end_header'), break; end
end

all_fields = strtrim(strsplit(fields_line, ','));
col_depth = find(strcmp(all_fields, 'depth'));
col_date  = find(strcmp(all_fields, 'date'));
dvsd_mask = strncmp(all_fields, 'PSD_DVSD_', 9);
col_dvsd  = find(dvsd_mask);
n_uvp     = numel(col_dvsd);

n_cols = numel(all_fields);
fmt    = repmat('%s ', 1, n_cols);
C      = textscan(fid, fmt, 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);
raw = C{1};

% parse depth, date, DVSD for surface rows only
depth_all = cellfun(@str2double, raw(:, col_depth));
date_all  = cellfun(@str2double, raw(:, col_date));
surf_mask = depth_all <= 5;

depth_surf = depth_all(surf_mask);
date_surf  = date_all(surf_mask);
raw_surf   = raw(surf_mask, :);
n_surf     = sum(surf_mask);

DVSD_surf = zeros(n_surf, n_uvp);
for j = 1:n_uvp
    vals = cellfun(@str2double, raw_surf(:, col_dvsd(j)));
    vals(abs(vals - missing_val) < 1) = NaN;
    DVSD_surf(:, j) = vals;
end

% bin widths [um] (same as parse_uvp)
d_bounds = [50.8, 64, 80.6, 102, 128, 161, 203, 256, 323, 406, 512, 645, ...
            813, 1020, 1290, 1630, 2050, 2580, 3250, 4100, 5160, 6500, ...
            8190, 10300, 13000, 16400, 20600, 26000];
dw = diff(d_bounds);

% phi per row [cm^3/cm^3]
phi_surf = DVSD_surf .* dw * 1e-9;

% --- average per date ---
unique_dates = unique(date_surf);
n_obs = numel(unique_dates);
phi_by_date = NaN(n_obs, n_uvp);
for i = 1:n_obs
    rows = date_surf == unique_dates(i);
    phi_by_date(i, :) = mean_no_nan(phi_surf(rows, :), 1);
end

% --- map to model bins ---
grid   = cfg.derive();
n_sec  = cfg.n_sections;
r_cm   = (0.75 / pi * grid.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;

% extract UVP bin centers from field names
d_uvp_um = zeros(1, n_uvp);
for i = 1:n_uvp
    tok = regexp(all_fields{col_dvsd(i)}, 'PSD_DVSD_(\d+)umsize', 'tokens');
    d_uvp_um(i) = str2double(tok{1}{1});
end

% nearest model bin for each UVP bin
bin_map = zeros(1, n_uvp);
for i = 1:n_uvp
    [~, bin_map(i)] = min(abs(d_model_um - d_uvp_um(i)));
end

phi_model_by_date = zeros(n_obs, n_sec);
for i = 1:n_uvp
    k = bin_map(i);
    vals = phi_by_date(:, i);
    vals(isnan(vals)) = 0;
    phi_model_by_date(:, k) = phi_model_by_date(:, k) + vals;
end

% power-law fill for model bins below 100 um (no UVP data there)
% fit log-log line to UVP 100-400 um range, extrapolate down
mask_small = d_model_um < 100;
fit_range  = d_uvp_um >= 100 & d_uvp_um <= 400;
for id = 1:n_obs
    phi_row = phi_by_date(id, :);
    phi_row(isnan(phi_row)) = 0;
    ok = fit_range & phi_row > 0;
    if sum(ok) >= 2
        p = polyfit(log10(d_uvp_um(ok)), log10(phi_row(ok)), 1);
        phi_fill = 10 .^ polyval(p, log10(d_model_um(mask_small)'));
        phi_fill(~isfinite(phi_fill) | phi_fill < 0) = 0;
        phi_model_by_date(id, mask_small) = phi_fill;
    end
end

% --- build daily time series (fill gaps) ---
% convert YYYYMMDD to day number from first date
first_date = unique_dates(1);
last_date  = unique_dates(end);

% convert to MATLAB datenum for arithmetic
dn_obs  = datenum(num2str(unique_dates), 'yyyymmdd');
dn_all  = (dn_obs(1) : dn_obs(end))';
n_days  = numel(dn_all);
dates_all = str2double(cellstr(datestr(dn_all, 'yyyymmdd')));

% interpolate each model bin
phi_daily = zeros(n_days, n_sec);
has_data  = false(n_days, 1);

% mark which days have real data
for i = 1:n_obs
    idx = find(dn_all == dn_obs(i));
    has_data(idx) = true;
end

% interpolate (pchip preserves positivity better than linear for sparse data)
for k = 1:n_sec
    y_obs = phi_model_by_date(:, k);
    if all(y_obs == 0)
        continue;
    end
    phi_daily(:, k) = max(0, interp1(dn_obs, y_obs, dn_all, 'pchip', 'extrap'));
end

% --- pack output ---
daily.day_num  = (1:n_days)';
daily.dates    = dates_all;
daily.phi      = phi_daily;       % n_days x n_sec
daily.has_data = has_data;
daily.d_model  = d_model_um(:)';
daily.n_days   = n_days;

end

function y = mean_no_nan(x, dim)
good = ~isnan(x);
x(~good) = 0;
n = sum(good, dim);
y = sum(x, dim) ./ max(n, 1);
y(n == 0) = NaN;
end

function uvp = parse_uvp(sb_file)
% PARSE_UVP  Read a UVP differential SeaBASS .sb file.
%
% Usage:
%   uvp = parse_uvp(sb_file)
%
% Output struct fields:
%   uvp.d_um      - UVP bin centers [um], 1 x 27
%   uvp.d_bounds  - bin boundaries [um], 1 x 28
%   uvp.dw        - bin widths [um], 1 x 27
%   uvp.depth_m   - unique depth levels [m], n_d x 1
%   uvp.DNSD      - mean DNSD over all casts [#/m^3/um], n_d x 27
%   uvp.DVSD      - mean DVSD over all casts [uL/m^3/um], n_d x 27
%   uvp.N         - number concentration per bin [#/m^3], n_d x 27
%   uvp.phi       - biovolume concentration [cm^3/cm^3], n_d x 27
%   uvp.date      - dates of each raw row (YYYYMMDD)
%   uvp.n_casts   - number of unique casts used
%
% Notes:
%   - Missing values (-9999) replaced with NaN before averaging.
%   - DNSD averaged over all casts at each depth bin.
%   - phi uses UVP DVSD directly:
%       DVSD [uL/m^3/um] * dw [um] * 1e-9 = cm^3/cm^3

% --- read header ---
fid = fopen(sb_file, 'r');
if fid < 0
    error('Cannot open file: %s', sb_file);
end

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
    if strcmp(strtrim(line), '/end_header')
        break;
    end
end

% --- parse field names ---
all_fields = strsplit(fields_line, ',');
all_fields = strtrim(all_fields);

% find column indices for depth and DNSD bins
col_depth = find(strcmp(all_fields, 'depth'));
col_date  = find(strcmp(all_fields, 'date'));
dnsd_mask = strncmp(all_fields, 'PSD_DNSD_', 9);
vold_mask = strncmp(all_fields, 'PSD_DVSD_', 9);
col_dnsd  = find(dnsd_mask);
col_dvsd  = find(vold_mask);
n_bins    = numel(col_dnsd);

if n_bins == 0 || numel(col_dvsd) ~= n_bins
    error('Could not find matching DNSD and DVSD UVP columns.');
end

% extract size in um from field names like PSD_DNSD_57umsize
d_um = zeros(1, n_bins);
for i = 1:n_bins
    tok = regexp(all_fields{col_dnsd(i)}, 'PSD_DNSD_(\d+)umsize', 'tokens');
    d_um(i) = str2double(tok{1}{1});
end

% bin boundaries from header values, 28 values for 27 bins
d_bounds = [50.8, 64, 80.6, 102, 128, 161, 203, 256, 323, 406, 512, 645, ...
            813, 1020, 1290, 1630, 2050, 2580, 3250, 4100, 5160, 6500, ...
            8190, 10300, 13000, 16400, 20600, 26000];
dw = diff(d_bounds);   % bin widths [um]

% --- read data rows ---
n_cols = numel(all_fields);
fmt    = repmat('%s ', 1, n_cols);
C      = textscan(fid, fmt, 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);

raw = C{1};
n_rows = size(raw, 1);

% parse depth and date
depth_all = cellfun(@str2double, raw(:, col_depth));
date_all  = cellfun(@str2double, raw(:, col_date));

% parse DNSD matrix
DNSD_raw = zeros(n_rows, n_bins);
DVSD_raw = zeros(n_rows, n_bins);
for j = 1:n_bins
    vals = cellfun(@str2double, raw(:, col_dnsd(j)));
    vals(abs(vals - missing_val) < 1) = NaN;   % replace missing
    DNSD_raw(:, j) = vals;

    vals = cellfun(@str2double, raw(:, col_dvsd(j)));
    vals(abs(vals - missing_val) < 1) = NaN;
    DVSD_raw(:, j) = vals;
end

% --- average over casts at each depth ---
depth_levels = unique(depth_all);
n_d = numel(depth_levels);

DNSD_mean = NaN(n_d, n_bins);
DVSD_mean = NaN(n_d, n_bins);
for i = 1:n_d
    rows = depth_all == depth_levels(i);
    DNSD_mean(i, :) = mean_no_nan(DNSD_raw(rows, :), 1);
    DVSD_mean(i, :) = mean_no_nan(DVSD_raw(rows, :), 1);
end

% number concentration [#/m^3] = DNSD * bin_width
N = DNSD_mean .* dw;   % n_d x n_bins

% biovolume concentration [cm^3/cm^3]
% DVSD [uL/m^3/um] * dw [um] = uL/m^3; 1 uL/m^3 = 1e-9 cm^3/cm^3
phi = DVSD_mean .* dw * 1e-9;

% count unique casts (unique station+date combinations)
station_col = find(strcmp(all_fields, 'station'));
if ~isempty(station_col)
    stations = raw(:, station_col);
    cast_id = cell(size(stations));
    for i = 1:numel(stations)
        cast_id{i} = sprintf('%s_%08d', stations{i}, date_all(i));
    end
    n_casts = numel(unique(cast_id));
else
    n_casts = numel(unique(date_all));
end

% --- pack output ---
uvp.d_um     = d_um;
uvp.d_bounds = d_bounds;
uvp.dw       = dw;
uvp.depth_m  = depth_levels;
uvp.DNSD     = DNSD_mean;
uvp.DVSD     = DVSD_mean;
uvp.N        = N;
uvp.phi      = phi;
uvp.date     = date_all;
uvp.n_casts  = n_casts;

end

function y = mean_no_nan(x, dim)
% Mean without needing special toolboxes.
good = ~isnan(x);
x(~good) = 0;
n = sum(good, dim);
y = sum(x, dim) ./ max(n, 1);
y(n == 0) = NaN;
end

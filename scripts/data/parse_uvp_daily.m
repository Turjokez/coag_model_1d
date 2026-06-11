function uvpd = parse_uvp_daily(sb_file)
% PARSE_UVP_DAILY  Read UVP .sb file, grouped by cast date and depth.
%
% Usage:
%   uvpd = parse_uvp_daily(sb_file)
%
% Output struct:
%   uvpd.dates    - unique cast dates [YYYYMMDD], n_dates x 1
%   uvpd.depth_m  - unique depth levels [m], n_depths x 1
%   uvpd.phi      - biovolume [cm^3/cm^3], n_dates x n_depths x n_bins
%                   NaN where no data for that date/depth pair
%   uvpd.d_um     - UVP bin centers [um], 1 x n_bins
%   uvpd.d_bounds - bin boundaries [um], 1 x 28
%   uvpd.dw       - bin widths [um], 1 x n_bins
%
% Notes:
%   - Multiple rows with same date+depth are averaged.
%   - phi = DVSD [uL/m^3/um] * dw [um] * 1e-9  [cm^3/cm^3]

% --- read header ---
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
col_depth  = find(strcmp(all_fields, 'depth'));
col_date   = find(strcmp(all_fields, 'date'));
dvsd_mask  = strncmp(all_fields, 'PSD_DVSD_', 9);
col_dvsd   = find(dvsd_mask);
n_bins     = numel(col_dvsd);

if n_bins == 0
    error('No PSD_DVSD_* columns found.');
end

% read all rows
n_cols = numel(all_fields);
fmt    = repmat('%s ', 1, n_cols);
C      = textscan(fid, fmt, 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);
raw = C{1};
n_rows = size(raw, 1);

% parse depth, date, DVSD
depth_all = cellfun(@str2double, raw(:, col_depth));
date_all  = cellfun(@str2double, raw(:, col_date));

DVSD_raw = zeros(n_rows, n_bins);
for j = 1:n_bins
    vals = cellfun(@str2double, raw(:, col_dvsd(j)));
    vals(abs(vals - missing_val) < 1) = NaN;
    DVSD_raw(:, j) = vals;
end

% bin geometry (same as parse_uvp)
d_bounds = [50.8, 64, 80.6, 102, 128, 161, 203, 256, 323, 406, 512, 645, ...
            813, 1020, 1290, 1630, 2050, 2580, 3250, 4100, 5160, 6500, ...
            8190, 10300, 13000, 16400, 20600, 26000];
dw = diff(d_bounds);

% bin centers from field names
d_um = zeros(1, n_bins);
for i = 1:n_bins
    tok = regexp(all_fields{col_dvsd(i)}, 'PSD_DVSD_(\d+)umsize', 'tokens');
    d_um(i) = str2double(tok{1}{1});
end

% --- group by date x depth ---
unique_dates  = unique(date_all);
unique_depths = unique(depth_all);
n_dates  = numel(unique_dates);
n_depths = numel(unique_depths);

phi_out = NaN(n_dates, n_depths, n_bins);

for id = 1:n_dates
    for iz = 1:n_depths
        rows = date_all == unique_dates(id) & depth_all == unique_depths(iz);
        if ~any(rows), continue; end
        DVSD_slice = DVSD_raw(rows, :);   % average if multiple rows
        DVSD_mean  = mean_no_nan(DVSD_slice, 1);
        phi_out(id, iz, :) = DVSD_mean .* dw * 1e-9;
    end
end

% --- pack output ---
uvpd.dates    = unique_dates;
uvpd.depth_m  = unique_depths;
uvpd.phi      = phi_out;    % n_dates x n_depths x n_bins
uvpd.d_um     = d_um;
uvpd.d_bounds = d_bounds;
uvpd.dw       = dw;

end

function y = mean_no_nan(x, dim)
good = ~isnan(x);
x(~good) = 0;
n = sum(good, dim);
y = sum(x, dim) ./ max(n, 1);
y(n == 0) = NaN;
end

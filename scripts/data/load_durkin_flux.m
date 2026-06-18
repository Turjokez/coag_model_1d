function d = load_durkin_flux(filepath)
% load_durkin_flux  Load Durkin gel trap particle fluxes from SeaBASS file.
%
% Returns struct d:
%   d.depths   [n_dep x 1]      unique trap depths [m]
%   d.d_um     [n_bin x 1]      bin center diameters [um]
%   d.flux_agg [n_dep x n_bin]  aggregate (ID2) flux [particles m-2 d-1]
%   d.flux_fp  [n_dep x n_bin]  fecal pellet (ID3+ID6+ID7) flux
%
% Flux averaged over all available dates at each depth.
% Missing values (-9999) treated as 0.
%
% Column indices (from SeaBASS /fields header):
%   6  = depth [m]
%   12 = bin_diameter_center [um]
%   16 = flux_particles_2id  (aggregate)
%   20 = flux_particles_3id  (long fecal pellet)
%   24 = flux_particles_7id  (mini pellet)
%   28 = flux_particles_6id  (short fecal pellet)

fid = fopen(filepath, 'r');
if fid < 0
    error('Cannot open: %s', filepath);
end

% skip header block
while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && contains(line, '/end_header')
        break;
    end
end

% read data lines
raw = struct('depth', {}, 'date', {}, 'd_um', {}, ...
             'agg', {}, 'fp3', {}, 'fp7', {}, 'fp6', {});
idx = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line) || isempty(line) || line(1) == '!'
        continue;
    end
    p = strsplit(line, ',');
    if numel(p) < 28
        continue;
    end
    depth = str2double(p{6});
    date  = str2double(p{8});
    d_um  = str2double(p{12});
    agg   = str2double(p{16});
    fp3   = str2double(p{20});
    fp7   = str2double(p{24});
    fp6   = str2double(p{28});
    if isnan(depth) || isnan(d_um)
        continue;
    end
    % replace missing flag with 0
    agg(agg < 0) = 0;
    fp3(fp3 < 0) = 0;
    fp7(fp7 < 0) = 0;
    fp6(fp6 < 0) = 0;
    idx = idx + 1;
    raw(idx).depth = depth;
    raw(idx).date  = date;
    raw(idx).d_um  = d_um;
    raw(idx).agg   = agg;
    raw(idx).fp3   = fp3;
    raw(idx).fp7   = fp7;
    raw(idx).fp6   = fp6;
end
fclose(fid);

if idx == 0
    error('No data rows found in %s', filepath);
end

% unique depths and bin centers
all_depths = [raw.depth]';
all_d_um   = [raw.d_um]';
depths     = unique(round(all_depths));     % round to nearest m
d_um       = unique(all_d_um);

n_dep = numel(depths);
n_bin = numel(d_um);

% accumulate flux sums and counts per (depth, d_um)
sum_agg = zeros(n_dep, n_bin);
sum_fp  = zeros(n_dep, n_bin);
cnt     = zeros(n_dep, n_bin);

for i = 1:idx
    id = find(depths == round(raw(i).depth), 1);
    ib = find(d_um   == raw(i).d_um, 1);
    if isempty(id) || isempty(ib)
        continue;
    end
    sum_agg(id, ib) = sum_agg(id, ib) + raw(i).agg;
    sum_fp(id, ib)  = sum_fp(id, ib)  + raw(i).fp3 + raw(i).fp7 + raw(i).fp6;
    cnt(id, ib)     = cnt(id, ib) + 1;
end

cnt(cnt == 0) = NaN;
d.depths   = depths;
d.d_um     = d_um;
d.flux_agg = sum_agg ./ cnt;
d.flux_fp  = sum_fp  ./ cnt;

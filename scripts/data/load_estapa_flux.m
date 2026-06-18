function d = load_estapa_flux(stt_file, nbst_file)
% load_estapa_flux  Load Estapa STT + NBST POC flux from SeaBASS files.
%
% Returns struct d:
%   d.depths    [n_dep x 1]  unique trap depths [m]
%   d.flux_POC  [n_dep x 1]  mean POC flux [mg C m-2 d-1]
%   d.flux_POC_sd  [n_dep x 1]  std of replicates [mg C m-2 d-1]
%   d.flux_Th   [n_dep x 1]  mean Th-234 flux [dpm m-2 d-1]
%   d.source    [n_dep x 1]  cell array: 'STT' or 'NBST'
%
% Column indices (same for STT and NBST):
%   3  = depth [m]
%   30 = flux_POC  [mmol C m-2 d-1]  -> convert to mg C by *12
%   38 = flux_Th_234 [dpm m-2 d-1]
%
% Missing values (-9999) excluded.
% Averaged over all replicates (A/B/C) and epochs at each depth.

rows_stt  = parse_file(stt_file,  'STT');
rows_nbst = parse_file(nbst_file, 'NBST');
rows = [rows_stt(:); rows_nbst(:)];

% unique depths
all_z   = [rows.depth]';
all_poc = [rows.poc]';   % mmol C m-2 d-1
all_th  = [rows.th]';    % dpm m-2 d-1
all_src = {rows.source}';

depths  = unique(all_z);
n_dep   = numel(depths);

d.depths       = depths;
d.flux_POC     = zeros(n_dep, 1);
d.flux_POC_sd  = zeros(n_dep, 1);
d.flux_Th      = zeros(n_dep, 1);
d.source       = cell(n_dep, 1);

for i = 1:n_dep
    mask = all_z == depths(i);
    poc  = all_poc(mask);
    th   = all_th(mask);
    src  = all_src(mask);
    d.flux_POC(i)    = mean(poc,  'omitnan') * 12;   % -> mg C m-2 d-1
    d.flux_POC_sd(i) = std(poc,   'omitnan') * 12;
    d.flux_Th(i)     = mean(th,   'omitnan');
    d.source{i}      = src{1};
end

% ---------------------------------------------------------------
function rows = parse_file(filepath, src_label)
fid = fopen(filepath, 'r');
if fid < 0, error('Cannot open: %s', filepath); end

while ~feof(fid)
    line = fgetl(fid);
    if ischar(line) && contains(line, '/end_header'), break; end
end

rows = struct('depth', {}, 'poc', {}, 'th', {}, 'source', {});
idx  = 0;
while ~feof(fid)
    line = fgetl(fid);
    if ~ischar(line) || isempty(strtrim(line)) || line(1) == '!'
        continue;
    end
    p = strsplit(line, ',');
    if numel(p) < 38, continue; end
    depth = str2double(p{3});
    poc   = str2double(p{30});
    th    = str2double(p{38});
    if isnan(depth), continue; end
    if poc <= -9999, poc = NaN; end
    if th  <= -9999, th  = NaN; end
    if isnan(poc) && isnan(th), continue; end
    idx = idx + 1;
    rows(idx).depth  = depth;
    rows(idx).poc    = poc;
    rows(idx).th     = th;
    rows(idx).source = src_label;
end
fclose(fid);

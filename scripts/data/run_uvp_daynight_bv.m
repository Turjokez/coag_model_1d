% run_uvp_daynight_bv.m
%
% DVM test: compare mean UVP biovolume BV(z) for all day casts
% vs all night casts.
%
% DVM prediction: if zooplankton migrate to depth during the day
% and produce fecal pellets there, daytime BV at 350-500 m should
% exceed nighttime BV.
%
% Day   = UTC 06:00-20:00  (sun up at ~49N in May)
% Night = UTC 20:00-06:00
% Bins: 100-2000 um only (UVP reliable range).

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);

uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Read UVP file
% ---------------------------------------------------------------
fid = fopen(uvp_file, 'r');
fields_line = '';
while true
    line = fgetl(fid);
    if ~ischar(line), break; end
    if startsWith(line, '/fields='), fields_line = extractAfter(line, '/fields='); end
    if strcmp(strtrim(line), '/end_header'), break; end
end
all_fields = strtrim(strsplit(fields_line, ','));
col_depth = find(strcmp(all_fields, 'depth'));
col_time  = find(strcmp(all_fields, 'time'));
dvsd_mask = strncmp(all_fields, 'PSD_DVSD_', 9);
col_dvsd  = find(dvsd_mask);
n_uvp_all = numel(col_dvsd);

fmt = repmat('%s ', 1, numel(all_fields));
C   = textscan(fid, fmt, 'Delimiter', ',', 'CollectOutput', true);
fclose(fid);
raw = C{1};

depth_all = cellfun(@str2double, raw(:, col_depth));
hour_all  = cellfun(@(t) str2double(t(1:2)), raw(:, col_time));

DVSD = zeros(size(raw, 1), n_uvp_all);
for j = 1:n_uvp_all
    vals = cellfun(@str2double, raw(:, col_dvsd(j)));
    vals(vals < -999) = NaN;
    DVSD(:, j) = vals;
end

% ---------------------------------------------------------------
% 2. Bin geometry — filter 100-2000 um
% ---------------------------------------------------------------
d_bounds = [50.8, 64, 80.6, 102, 128, 161, 203, 256, 323, 406, 512, 645, ...
            813, 1020, 1290, 1630, 2050, 2580, 3250, 4100, 5160, 6500, ...
            8190, 10300, 13000, 16400, 20600, 26000];
d_centers = sqrt(d_bounds(1:end-1) .* d_bounds(2:end));
dw_um     = diff(d_bounds);
mask_uvp  = d_centers >= 100 & d_centers < 2000;
dw_filt   = dw_um(mask_uvp);        % bin widths in um

% particle volume per bin [cm^3], assume sphere
r_uvp_cm  = (d_centers(mask_uvp) / 2) / 1e4;
av_vol    = (4/3) * pi * r_uvp_cm.^3;   % cm^3

% ---------------------------------------------------------------
% 3. Depth grid
% ---------------------------------------------------------------
z_edges = 50:25:500;
z_mid   = (z_edges(1:end-1) + z_edges(2:end)) / 2;
n_z     = numel(z_mid);
n_bins  = sum(mask_uvp);

% ---------------------------------------------------------------
% 4. Separate rows into day (06-20 UTC) and night (20-06 UTC)
% ---------------------------------------------------------------
is_day   = hour_all >= 6 & hour_all < 20;
is_night = ~is_day;

fprintf('Day rows: %d   Night rows: %d\n', sum(is_day), sum(is_night));

% ---------------------------------------------------------------
% 5. Compute mean BV(z) and N(z) for day and night
% ---------------------------------------------------------------
BV_day   = NaN(n_z, 1);
BV_night = NaN(n_z, 1);
N_day    = NaN(n_z, 1);
N_night  = NaN(n_z, 1);

for iz = 1:n_z
    in_layer = depth_all >= z_edges(iz) & depth_all < z_edges(iz+1);

    % day
    rows_d = in_layer & is_day;
    if any(rows_d)
        phi_d = DVSD(rows_d, mask_uvp);   % [nrows x nbins], ppmV/um
        phi_d(phi_d < 0) = NaN;
        % convert DVSD [ppmV/um] to phi [ppmV per bin]: multiply by dw_um
        phi_ppmv = bsxfun(@times, phi_d, dw_filt);
        BV_day(iz) = mean(mean(phi_ppmv, 2, 'omitnan'), 'omitnan');
        % N = phi_ppmv / av_vol [cm^3], then sum over bins
        N_day(iz)  = mean(sum(bsxfun(@rdivide, phi_ppmv * 1e-6, av_vol), 2, 'omitnan'), 'omitnan');
    end

    % night
    rows_n = in_layer & is_night;
    if any(rows_n)
        phi_n = DVSD(rows_n, mask_uvp);
        phi_n(phi_n < 0) = NaN;
        phi_ppmv = bsxfun(@times, phi_n, dw_filt);
        BV_night(iz) = mean(mean(phi_ppmv, 2, 'omitnan'), 'omitnan');
        N_night(iz)  = mean(sum(bsxfun(@rdivide, phi_ppmv * 1e-6, av_vol), 2, 'omitnan'), 'omitnan');
    end
end

% ---------------------------------------------------------------
% 6. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 14 9], 'Color', 'white');

ax1 = subplot(1, 2, 1);
plot(BV_day,   z_mid, 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
plot(BV_night, z_mid, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2);
set(ax1, 'YDir', 'reverse', 'YLim', [50 500], 'XScale', 'log');
xlabel('Biovolume (ppmV)');
ylabel('Depth (m)');
legend('Day (06-20 UTC)', 'Night (20-06 UTC)', 'Location', 'southeast', 'FontSize', 7);
title('a) BV', 'FontWeight', 'normal');

ax2 = subplot(1, 2, 2);
plot(N_day,   z_mid, 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
plot(N_night, z_mid, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2);
set(ax2, 'YDir', 'reverse', 'YLim', [50 500], 'XScale', 'log');
xlabel('Particle number (cm^{-3})');
legend('Day', 'Night', 'Location', 'southeast', 'FontSize', 7);
title('b) N', 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'uvp_daynight_bv.png'));
fprintf('Saved uvp_daynight_bv.png\n');

% ---------------------------------------------------------------
% 7. Print ratio table
% ---------------------------------------------------------------
fprintf('\nDepth   BV_day/BV_night   N_day/N_night\n');
for iz = 1:n_z
    fprintf('%5.0f m   %6.2f            %6.2f\n', ...
        z_mid(iz), BV_day(iz)/BV_night(iz), N_day(iz)/N_night(iz));
end

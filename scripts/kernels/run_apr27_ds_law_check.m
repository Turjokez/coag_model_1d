% run_apr27_ds_law_check
% Short note:
% 1. check if beta_ds uses the selected sinking law
% 2. compare one fixed-size curve and one size-size map
% 3. keep this inside the 1-D testing folder

clear;
clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
tab_dir = fullfile(repo_root, 'output', 'tables');
log_dir = fullfile(repo_root, 'output', 'logs');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(tab_dir, 'dir')
    mkdir(tab_dir);
end
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end

law_names = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
law_labels = {'current', 'kriest 8', 'kriest 9', 'siegel 2025'};
mode_tag = 'direct_sinking_law';
d_um = logspace(0, 4, 320);
d_cm = d_um .* 1e-4;
d1_um = 1.0;
d1_cm = d1_um .* 1e-4;

fig1 = figure('Color', 'w', 'Position', [120 120 720 520]);
ax1 = axes(fig1);
hold(ax1, 'on');
cols = lines(numel(law_names));
for i = 1:numel(law_names)
    beta = local_beta_diff_sed_from_law(d1_cm, d_cm, law_names{i});
    plot(ax1, d_um, beta, 'LineWidth', 1.5, 'Color', cols(i, :), ...
        'DisplayName', law_labels{i});
end
set(ax1, 'XScale', 'log', 'YScale', 'log');
xlabel(ax1, 'Partner size (um)');
ylabel(ax1, 'beta ds (cm^3 s^{-1})');
title(ax1, ['Fixed 1 um particle (' mode_tag ')']);
legend(ax1, 'Location', 'northwest', 'Box', 'off');
grid(ax1, 'on');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
local_save_figure(fig1, fullfile(fig_dir, 'apr27_ds_law_fixed_1um.png'));
close(fig1);

[D1, D2] = ndgrid(d_cm, d_cm);
map_vals = cell(numel(law_names), 1);
all_log_vals = [];
for i = 1:numel(law_names)
    beta = local_beta_diff_sed_from_law(D1, D2, law_names{i});
    beta(beta <= 0) = NaN;
    map_vals{i} = log10(beta);
    all_log_vals = [all_log_vals; map_vals{i}(isfinite(map_vals{i}))]; %#ok<AGROW>
end
clim = [min(all_log_vals), max(all_log_vals)];

fig2 = figure('Color', 'w', 'Position', [120 120 950 760]);
tl = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(law_names)
    ax = nexttile(tl, i);
    imagesc(ax, log10(d_um), log10(d_um), map_vals{i});
    axis(ax, 'xy');
    caxis(ax, clim);
    xlabel(ax, 'log10 d1 (um)');
    ylabel(ax, 'log10 d2 (um)');
    title(ax, law_labels{i});
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
end
colormap(fig2, parula);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'log10 beta ds';
title(tl, ['Differential-settling beta from selected law (' mode_tag ')']);
local_save_figure(fig2, fullfile(fig_dir, 'apr27_ds_law_beta_maps.png'));
close(fig2);

fig3 = make_powerlaw_normalization_check(d1_um, d_um, mode_tag);
local_save_figure(fig3, fullfile(fig_dir, 'apr27_ds_powerlaw_normalization_check.png'));
close(fig3);

sample_d2_um = [10; 100; 1000; 10000];
rows = struct('law', {}, 'd1_um', {}, 'd2_um', {}, ...
    'w1_cm_s', {}, 'w2_cm_s', {}, 'dw_cm_s', {}, 'beta_cm3_s', {});
for i = 1:numel(law_names)
    for j = 1:numel(sample_d2_um)
        d2_cm = sample_d2_um(j) .* 1e-4;
        [beta, w1, w2] = local_beta_diff_sed_from_law(d1_cm, d2_cm, law_names{i});
        row = struct();
        row.law = string(law_names{i});
        row.d1_um = d1_um;
        row.d2_um = sample_d2_um(j);
        row.w1_cm_s = w1;
        row.w2_cm_s = w2;
        row.dw_cm_s = abs(w2 - w1);
        row.beta_cm3_s = beta;
        rows(end + 1) = row; %#ok<AGROW>
    end
end
T = struct2table(rows);
csv_path = fullfile(tab_dir, 'apr27_ds_law_sample_checks.csv');
writetable(T, csv_path);

log_path = fullfile(log_dir, 'apr27_ds_law_check.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'April 27 differential-settling law check\n\n');
fprintf(fid, 'Purpose:\n');
fprintf(fid, '- make the sinking-law choice explicit in beta_ds\n');
fprintf(fid, '- check a fixed 1 um particle against partner sizes\n');
fprintf(fid, '- check full size-size beta maps with the same color scale\n\n');
fprintf(fid, 'Kernel mode used for these figures: %s\n\n', mode_tag);
fprintf(fid, 'Files:\n');
fprintf(fid, '- output/figures/apr27_ds_law_fixed_1um.png\n');
fprintf(fid, '- output/figures/apr27_ds_law_beta_maps.png\n');
fprintf(fid, '- output/figures/apr27_ds_powerlaw_normalization_check.png\n');
fprintf(fid, '- output/tables/apr27_ds_law_sample_checks.csv\n\n');
fprintf(fid, 'Reading:\n');
fprintf(fid, '- beta_ds is calculated from beta = (pi/4) * (d1+d2)^2 * abs(w1-w2)\n');
fprintf(fid, '- w1 and w2 come from the selected named sinking law\n');
fprintf(fid, '- the power-law check shows why normalization can change which curve looks larger\n');
fclose(fid);

disp('Saved differential-settling law check:');
disp(fullfile(fig_dir, 'apr27_ds_law_fixed_1um.png'));
disp(fullfile(fig_dir, 'apr27_ds_law_beta_maps.png'));
disp(fullfile(fig_dir, 'apr27_ds_powerlaw_normalization_check.png'));
disp(csv_path);
disp(log_path);

function fig = make_powerlaw_normalization_check(d1_um, d_um, mode_tag)
exponents = [0.6, 1.0, 2.0];
cols = lines(numel(exponents));
d1_cm = d1_um .* 1e-4;
d_cm = d_um .* 1e-4;
dmax_um = max(d_um);

fig = figure('Color', 'w', 'Position', [120 120 1080 460]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
hold(ax1, 'on');
for i = 1:numel(exponents)
    b = exponents(i);
    w1 = (d1_um ./ dmax_um) .^ b;
    w2 = (d_um ./ dmax_um) .^ b;
    beta = local_beta_diff_sed(d1_cm, d_cm, w1, w2);
    plot(ax1, d_um, beta, 'LineWidth', 1.4, 'Color', cols(i, :), ...
        'DisplayName', sprintf('b = %.1f', b));
end
set(ax1, 'XScale', 'log', 'YScale', 'log');
xlabel(ax1, 'Partner size (um)');
ylabel(ax1, 'beta ds (arb.)');
title(ax1, 'All speeds equal at 10000 um');
legend(ax1, 'Location', 'northwest', 'Box', 'off');
grid(ax1, 'on');

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
for i = 1:numel(exponents)
    b = exponents(i);
    w1 = d1_um .^ b;
    w2 = d_um .^ b;
    beta = local_beta_diff_sed(d1_cm, d_cm, w1, w2);
    plot(ax2, d_um, beta, 'LineWidth', 1.4, 'Color', cols(i, :), ...
        'DisplayName', sprintf('b = %.1f', b));
end
set(ax2, 'XScale', 'log', 'YScale', 'log');
xlabel(ax2, 'Partner size (um)');
ylabel(ax2, 'beta ds (arb.)');
title(ax2, 'Same prefactor for all b');
legend(ax2, 'Location', 'northwest', 'Box', 'off');
grid(ax2, 'on');

title(tl, ['Power-law beta check for fixed 1 um particle (' mode_tag ')']);
end

function [beta, w1, w2] = local_beta_diff_sed_from_law(d1_cm, d2_cm, law_name)
w1 = local_sinking_speed_named(d1_cm, law_name);
w2 = local_sinking_speed_named(d2_cm, law_name);
beta = local_beta_diff_sed(d1_cm, d2_cm, w1, w2);
end

function beta = local_beta_diff_sed(d1_cm, d2_cm, w1, w2)
beta = (pi/4) .* (d1_cm + d2_cm) .* (d1_cm + d2_cm) .* abs(w1 - w2);
end

function w_cm_s = local_sinking_speed_named(d_cm, law_name)
cfg = SimulationConfig();
law = lower(string(law_name));
switch law
    case "current"
        w_cm_s = local_current_law(d_cm, cfg);
    case "kriest_8"
        w_cm_s = (66 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "kriest_9"
        w_cm_s = (132 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "siegel_2025"
        d_mm = d_cm .* 10.0;
        w_cm_s = (20.2 .* (d_mm .^ 0.67) .* 100) ./ cfg.day_to_sec;
    otherwise
        error('Unknown law: %s', char(law_name));
end
w_cm_s(~isfinite(w_cm_s)) = 0;
w_cm_s(w_cm_s < 0) = 0;
end

function w_cm_s = local_current_law(d_cm, cfg)
r_v = 0.5 .* d_cm;
w_setcon = KernelLibrary.currentSetcon(cfg);
r_i = KernelLibrary.conservativeToFractalRadius(r_v, cfg);
w_cm_s = w_setcon .* (r_v .^ 3) ./ max(r_i, realmin);
end

function local_save_figure(fig_handle, fig_path)
[fig_parent, ~, ~] = fileparts(fig_path);
if ~exist(fig_parent, 'dir')
    mkdir(fig_parent);
end
set(fig_handle, 'PaperPositionMode', 'auto');
try
    exportgraphics(fig_handle, fig_path, 'Resolution', 220);
catch
    saveas(fig_handle, fig_path);
end
end

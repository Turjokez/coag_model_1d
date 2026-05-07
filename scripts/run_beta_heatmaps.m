% run_beta_heatmaps
% Steps:
% 1. build the size grid
% 2. compute beta on the grid
% 3. save simple heatmaps

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

diam_cm = logspace(-4, 0, 220);
law_names = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};

fig = plot_beta_heatmap(diam_cm, law_names);
out_path = fullfile(fig_dir, 'beta_heatmaps_named_laws.png');
save_figure(fig, out_path);

disp('Saved figure:');
disp(out_path);

function fig = plot_beta_heatmap(diam_cm, law_names)
diam_um = diam_cm .* 1e4;
[D1, D2] = ndgrid(diam_cm, diam_cm);

vals = cell(numel(law_names), 1);
all_vals = [];
for i = 1:numel(law_names)
    beta = local_beta_diff_sed_from_law(D1, D2, law_names{i});
    beta(beta <= 0) = NaN;
    vals{i} = log10(beta);
    all_vals = [all_vals; vals{i}(isfinite(vals{i}))]; %#ok<AGROW>
end

clim = [min(all_vals), max(all_vals)];
labels = {'current', 'kriest 8', 'kriest 9', 'siegel 2025'};

fig = figure('Color', 'w', 'Position', [100 100 900 730]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(law_names)
    ax = nexttile(tl, i);
    imagesc(ax, log10(diam_um), log10(diam_um), vals{i});
    axis(ax, 'xy');
    caxis(ax, clim);
    xlabel(ax, 'log10 d1 (um)');
    ylabel(ax, 'log10 d2 (um)');
    title(ax, labels{i});
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
end

colormap(fig, parula);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'log10 beta ds';
title(tl, 'Differential-settling beta from named laws');
end

function [beta, w1, w2] = local_beta_diff_sed_from_law(d1_cm, d2_cm, law_name)
w1 = local_sinking_speed_named(d1_cm, law_name);
w2 = local_sinking_speed_named(d2_cm, law_name);
beta = (pi / 4.0) .* (d1_cm + d2_cm) .* (d1_cm + d2_cm) .* abs(w1 - w2);
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
setcon = KernelLibrary.currentSetcon(cfg);
r_i = KernelLibrary.conservativeToFractalRadius(r_v, cfg);
w_cm_s = setcon .* (r_v .^ 3) ./ max(r_i, realmin);
end

function save_figure(fig_handle, fig_path)
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

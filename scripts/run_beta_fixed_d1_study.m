% run_beta_fixed_d1_study
% Steps:
% 1. choose one fixed size
% 2. sweep the second size
% 3. compare kernel curves

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Main settings. Edit a here if you want a different prefactor.
a = 1.0;
b_vals = 0.6:0.2:2.0;
d2_cm = logspace(-4, 0, 400);
d1_um_list = [1, 500, 1000];
d1_cm_list = d1_um_list * 1e-4;

for i = 1:numel(d1_cm_list)
    d1_cm = d1_cm_list(i);
    fig = plot_beta_fixed_d1(d1_cm, d2_cm, b_vals, a);
    out_name = sprintf('beta_fixed_d1_%s.png', make_d1_tag(d1_cm));
    out_path = fullfile(fig_dir, out_name);
    save_figure(fig, out_path);
    close(fig);
    disp('Saved figure:');
    disp(out_path);
end

function tag = make_d1_tag(d1_cm)
d1_um = d1_cm * 1e4;
if abs(d1_um - 1000) < 1e-12
    tag = '1mm';
else
    tag = sprintf('%gum', d1_um);
end
end

function fig = plot_beta_fixed_d1(d1_cm, d2_cm, b_vals, a)
d2_um = d2_cm .* 1e4;
cols = parula(numel(b_vals));

fig = figure('Color', 'w', 'Position', [100 100 760 500]);
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(b_vals)
    w1 = sinking_speed_powerlaw(d1_cm, a, b_vals(i));
    w2 = sinking_speed_powerlaw(d2_cm, a, b_vals(i));
    beta = local_beta_diff_sed(d1_cm, d2_cm, w1, w2);
    beta(beta <= 0) = NaN;
    plot(ax, d2_um, beta, 'LineWidth', 1.8, 'Color', cols(i, :), ...
        'DisplayName', sprintf('b = %.1f', b_vals(i)));
end

set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, 'Partner particle size, d2 (um)');
ylabel(ax, 'Differential-settling beta (comparison units)');
title(ax, sprintf('Fixed particle, d1 = %g um', d1_cm .* 1e4));
legend(ax, 'Location', 'northwest', 'NumColumns', 2, 'Box', 'off');
grid(ax, 'on');
ax.LineWidth = 1.0;
ax.FontSize = 11;
end

function beta = local_beta_diff_sed(d1_cm, d2_cm, w1, w2)
beta = (pi / 4.0) .* (d1_cm + d2_cm) .* (d1_cm + d2_cm) .* abs(w1 - w2);
end

function w = sinking_speed_powerlaw(d_cm, a, b)
w = a .* (d_cm .^ b);
w(~isfinite(w)) = 0;
w(w < 0) = 0;
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

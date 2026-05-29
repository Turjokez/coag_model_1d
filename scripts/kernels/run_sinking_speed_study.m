% run_sinking_speed_study
% Steps:
% 1. set paths
% 2. define the sweep
% 3. call the plotting helper

clear;
clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
tab_dir = fullfile(repo_root, 'output', 'tables');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(tab_dir, 'dir')
    mkdir(tab_dir);
end

% Main settings. Edit a here if you want a different prefactor.
a = 1.0;
b_vals = 0.6:0.2:2.0;
diam_cm = logspace(-4, 0, 400);

fig = plot_sinking_speed_sweep(diam_cm, b_vals, a);
fig_path = fullfile(fig_dir, 'sinking_speed_exponent_sweep.png');
save_figure(fig, fig_path);

sample_um = [1; 10; 100; 1000; 10000];
sample_cm = sample_um * 1e-4;
speed_mat = zeros(numel(sample_um), numel(b_vals));

for i = 1:numel(b_vals)
    speed_mat(:, i) = sinking_speed_powerlaw(sample_cm, a, b_vals(i));
end

var_names = [{'diameter_um', 'diameter_cm'}, make_b_names(b_vals)];
T = array2table([sample_um, sample_cm, speed_mat], 'VariableNames', var_names);

csv_path = fullfile(tab_dir, 'sinking_speed_exponent_sweep_summary.csv');
writetable(T, csv_path);

disp('Saved figure:');
disp(fig_path);
disp('Saved summary table:');
disp(csv_path);
disp('Representative sinking speeds:');
disp(T);

function names = make_b_names(b_vals)
names = cell(1, numel(b_vals));
for j = 1:numel(b_vals)
    txt = sprintf('b_%.1f', b_vals(j));
    txt = strrep(txt, '.', 'p');
    names{j} = matlab.lang.makeValidName(txt);
end
end

function fig = plot_sinking_speed_sweep(diam_cm, b_vals, a)
diam_um = diam_cm .* 1e4;
cols = parula(numel(b_vals));

fig = figure('Color', 'w', 'Position', [100 100 760 500]);
ax = axes(fig);
hold(ax, 'on');

for i = 1:numel(b_vals)
    w = sinking_speed_powerlaw(diam_cm, a, b_vals(i));
    plot(ax, diam_um, w, 'LineWidth', 1.8, 'Color', cols(i, :), ...
        'DisplayName', sprintf('b = %.1f', b_vals(i)));
end

set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, 'Particle diameter, d (um)');
ylabel(ax, 'Relative sinking speed, w');
legend(ax, 'Location', 'southeast', 'NumColumns', 2, 'Box', 'off');
grid(ax, 'on');
ax.LineWidth = 1.0;
ax.FontSize = 11;
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

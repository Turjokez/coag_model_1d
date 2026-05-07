% run_1d_step1_sinking_law_compare
% Steps:
% 1. compare named sinking laws in one simple setup
% 2. use the trusted advection-only upwind baseline
% 3. save simple figures and one short summary

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');
tab_dir = fullfile(repo_root, 'output', 'tables');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end
if ~exist(tab_dir, 'dir')
    mkdir(tab_dir);
end

law_names = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
law_labels = {'current', 'kriest 8', 'kriest 9', 'siegel 2025'};
plot_cols = lines(numel(law_names));

diam_plot_um = logspace(0, 4, 300);
diam_plot_cm = diam_plot_um .* 1e-4;
speed_plot_m_day = zeros(numel(diam_plot_um), numel(law_names));
travel_plot_day = zeros(numel(diam_plot_um), numel(law_names));

for i = 1:numel(law_names)
    speed_cm_s = sinking_speed_named(diam_plot_cm, law_names{i});
    speed_m_day = speed_cm_s .* 0.01 .* 86400.0;
    speed_plot_m_day(:, i) = speed_m_day;
    travel_plot_day(:, i) = 1000.0 ./ speed_m_day;
end

fig1 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax1 = axes(fig1);
hold(ax1, 'on');
for i = 1:numel(law_names)
    plot(ax1, diam_plot_um, speed_plot_m_day(:, i), 'LineWidth', 1.4, ...
        'Color', plot_cols(i, :), 'DisplayName', law_labels{i});
end
set(ax1, 'XScale', 'log', 'YScale', 'log');
xlabel(ax1, 'Diameter (um)');
ylabel(ax1, 'Sinking speed (m/day)');
title(ax1, 'Sinking speed to use in the 1-D model');
legend(ax1, 'Location', 'northwest', 'Box', 'off');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, 'step1_sinking_speed_laws.png'));
close(fig1);

fig2 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax2 = axes(fig2);
hold(ax2, 'on');
for i = 1:numel(law_names)
    plot(ax2, diam_plot_um, travel_plot_day(:, i), 'LineWidth', 1.4, ...
        'Color', plot_cols(i, :), 'DisplayName', law_labels{i});
end
set(ax2, 'XScale', 'log', 'YScale', 'log');
xlabel(ax2, 'Diameter (um)');
ylabel(ax2, 'Travel time to 1000 m (day)');
title(ax2, 'Travel time to 1000 m, t = 1000 / w');
legend(ax2, 'Location', 'southwest', 'Box', 'off');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, 'step1_travel_time_laws.png'));
close(fig2);

size_um = [100; 500; 1000; 3000];
size_cm = size_um .* 1e-4;
rows = struct('law', {}, 'size_um', {}, 'speed_m_day', {}, ...
    'expected_day', {}, 'peak_day', {}, 'error_pct', {});
bottom_curves_raw = cell(numel(law_names), 1);
t_day_raw = cell(numel(law_names), 1);

for i = 1:numel(law_names)
    speed_cm_s = sinking_speed_named(size_cm, law_names{i});
    speed_m_s = speed_cm_s .* 0.01;

    cfg = struct();
    cfg.z_max_m = 1000.0;
    cfg.dz_m = 5.0;
    cfg.dt_s = 0.90 .* cfg.dz_m ./ max(speed_m_s);
    cfg.t_max_s = 1.20 .* max(cfg.z_max_m ./ speed_m_s);
    cfg.size_um = size_um;
    cfg.speed_m_s = speed_m_s;
    cfg.pulse_amp = 1.0;
    cfg.law_name = law_names{i};
    cfg.scheme = 'upwind';

    sim = solve_advection_only(cfg);
    summary = local_travel_summary(sim);

    i_rep = find(size_um == 1000, 1, 'first');
    t_day_raw{i} = sim.t_s ./ 86400.0;
    bottom_curves_raw{i} = sim.bottom_signal(:, i_rep);

    for j = 1:height(summary)
        row = struct();
        row.law = string(law_names{i});
        row.size_um = summary.size_um(j);
        row.speed_m_day = summary.speed_m_day(j);
        row.expected_day = summary.expected_day(j);
        row.peak_day = summary.peak_day(j);
        row.error_pct = summary.error_pct(j);
        rows(end + 1) = row; %#ok<AGROW>
    end
end

t_max_day = 0;
for i = 1:numel(law_names)
    t_max_day = max(t_max_day, max(t_day_raw{i}));
end
t_day_ref = linspace(0, t_max_day, 2500).';
bottom_curves = zeros(numel(t_day_ref), numel(law_names));
for i = 1:numel(law_names)
    bottom_curves(:, i) = interp1(t_day_raw{i}, bottom_curves_raw{i}, ...
        t_day_ref, 'linear', 0.0);
end

fig3 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax3 = axes(fig3);
hold(ax3, 'on');
for i = 1:numel(law_names)
    plot(ax3, t_day_ref, bottom_curves(:, i), 'LineWidth', 1.4, ...
        'Color', plot_cols(i, :), 'DisplayName', law_labels{i});
end
xlabel(ax3, 'Time (day)');
ylabel(ax3, 'Particles at 1000 m');
title(ax3, 'Arrival at 1000 m for 1 mm');
legend(ax3, 'Location', 'best', 'Box', 'off');
ax3.LineWidth = 1.0;
ax3.FontSize = 11;
save_figure(fig3, fullfile(fig_dir, 'step1_bottom_signal_1mm_laws.png'));
close(fig3);

T = struct2table(rows);
csv_path = fullfile(tab_dir, 'step1_sinking_law_compare_summary.csv');
writetable(T, csv_path);

mean_err = groupsummary(T, 'law', 'mean', 'error_pct');
speed_1mm = T(T.size_um == 1000, {'law', 'speed_m_day', 'expected_day', 'peak_day', 'error_pct'});
log_path = fullfile(log_dir, 'step1_sinking_law_compare.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Step 1: named sinking-law compare in advection-only model\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- scheme = upwind\n');
fprintf(fid, '- depth = 1000 m\n');
fprintf(fid, '- sizes = 100, 500, 1000, 3000 um\n');
fprintf(fid, '- same pulse for each law\n\n');

fprintf(fid, 'Mean travel-time error by law:\n');
for i = 1:height(mean_err)
    fprintf(fid, '- %s: %.4f %%\n', char(mean_err.law(i)), mean_err.mean_error_pct(i));
end

fprintf(fid, '\n1 mm summary:\n');
for i = 1:height(speed_1mm)
    fprintf(fid, '- %s: speed = %.3f m/day, expected = %.3f day, peak = %.3f day, error = %.4f %%\n', ...
        char(speed_1mm.law(i)), speed_1mm.speed_m_day(i), speed_1mm.expected_day(i), ...
        speed_1mm.peak_day(i), speed_1mm.error_pct(i));
end
fclose(fid);

disp('Saved step 1 figures and summary:');
disp(fullfile(fig_dir, 'step1_sinking_speed_laws.png'));
disp(fullfile(fig_dir, 'step1_travel_time_laws.png'));
disp(fullfile(fig_dir, 'step1_bottom_signal_1mm_laws.png'));
disp(csv_path);
disp(log_path);

function out = local_travel_summary(sim)
t_day = sim.t_s(:) ./ 86400.0;
expected_day = sim.cfg.z_max_m ./ sim.speed_m_day;
ns = numel(sim.size_um);

peak_day = zeros(ns, 1);
for is = 1:ns
    y = sim.bottom_signal(:, is);
    [~, idx] = max(y);
    peak_day(is) = t_day(idx);
end

err_pct = 100.0 .* (peak_day - expected_day) ./ max(expected_day, realmin);

out = table(sim.size_um, sim.speed_m_day, expected_day, peak_day, err_pct, ...
    'VariableNames', {'size_um', 'speed_m_day', 'expected_day', 'peak_day', 'error_pct'});
end

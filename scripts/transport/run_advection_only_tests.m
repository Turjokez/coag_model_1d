% run_advection_only_tests
% Steps:
% 1. run both advection schemes
% 2. compare travel time and spreading
% 3. save a short summary

clear;
clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');
table_dir = fullfile(repo_root, 'output', 'tables');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end
if ~exist(table_dir, 'dir')
    mkdir(table_dir);
end

% Keep the first test simple.
% Use a few larger sizes so the travel time to 1000 m stays practical.
law_name = 'kriest_8';
size_um = [100; 500; 1000; 3000];
size_cm = size_um ./ 1e4;
speed_cm_s = local_named_speed(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

z_max_m = 1000.0;
dz_m = 5.0;
dt_s = 0.90 .* dz_m ./ max(speed_m_s);
travel_s = z_max_m ./ speed_m_s;
t_max_s = 1.20 .* max(travel_s);

cfg = struct();
cfg.z_max_m = z_max_m;
cfg.dz_m = dz_m;
cfg.dt_s = dt_s;
cfg.t_max_s = t_max_s;
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
cfg.pulse_amp = 1.0;
cfg.law_name = law_name;

schemes = {'upwind', 'lax_wendroff'};
res = repmat(struct('scheme', '', 'sim', [], 'val', [], 'csv_path', ''), 0, 1);

for i = 1:numel(schemes)
    cfg.scheme = schemes{i};
    sim = solve_advection_only(cfg);
    val = validate_travel_time(sim, fig_dir, log_dir);

    csv_path = fullfile(table_dir, sprintf('advection_only_%s_validation_summary.csv', schemes{i}));
    writetable(val.summary, csv_path);

    row = struct();
    row.scheme = schemes{i};
    row.sim = sim;
    row.val = val;
    row.csv_path = csv_path;
    res(end + 1, 1) = row; %#ok<AGROW>
end

save_comparison_figure(res, fig_dir);
summary_log = save_summary_log(res, log_dir);

disp('Saved comparison log:');
disp(summary_log);
for i = 1:numel(res)
    disp(res(i).csv_path);
end

function w = local_named_speed(diam_cm, law_name)
switch lower(string(law_name))
    case "current"
        w = sinking_speed_current(diam_cm);
    case "kriest_8"
        w = sinking_speed_kriest8(diam_cm);
    case "kriest_9"
        w = sinking_speed_kriest9(diam_cm);
    case "siegel_2025"
        w = sinking_speed_siegel2025(diam_cm);
    otherwise
        error('run_advection_only_tests:law', 'Unknown law: %s', law_name);
end
end

function save_comparison_figure(res, fig_dir)
fig = figure('Color', 'w', 'Position', [70 70 1100 860]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

rep_size_um = 1000;
i_rep = find(res(1).sim.size_um == rep_size_um, 1, 'first');
if isempty(i_rep)
    i_rep = min(2, numel(res(1).sim.size_um));
end

ax1 = nexttile(tl, 1);
hold(ax1, 'on');
for i = 1:numel(res)
    t_day = res(i).sim.t_s ./ 86400.0;
    y = res(i).sim.bottom_signal(:, i_rep);
    plot(ax1, t_day, y, 'LineWidth', 1.4, 'DisplayName', res(i).scheme);
end
xlabel(ax1, 'Time (day)');
ylabel(ax1, 'Signal at 1000 m');
title(ax1, sprintf('Bottom signal | %.1f mm', res(1).sim.size_um(i_rep) / 1000.0));
grid(ax1, 'on');
legend(ax1, 'Location', 'best');

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
for i = 1:numel(res)
    plot(ax2, res(i).val.summary.size_um, abs(res(i).val.summary.error_pct), '-o', ...
        'LineWidth', 1.3, 'MarkerSize', 6, 'DisplayName', res(i).scheme);
end
xlabel(ax2, 'Size (um)');
ylabel(ax2, '|Travel-time error| (%)');
title(ax2, 'Travel-time accuracy');
grid(ax2, 'on');
legend(ax2, 'Location', 'best');

ax3 = nexttile(tl, 3);
hold(ax3, 'on');
for i = 1:numel(res)
    plot(ax3, res(i).val.summary.size_um, res(i).val.summary.signal_width_day, '-o', ...
        'LineWidth', 1.3, 'MarkerSize', 6, 'DisplayName', res(i).scheme);
end
xlabel(ax3, 'Size (um)');
ylabel(ax3, 'Bottom-signal width (day)');
title(ax3, 'Artificial spreading');
grid(ax3, 'on');

ax4 = nexttile(tl, 4);
hold(ax4, 'on');
for i = 1:numel(res)
    plot(ax4, res(i).val.summary.size_um, res(i).val.summary.mass_drop_pct, '-o', ...
        'LineWidth', 1.3, 'MarkerSize', 6, 'DisplayName', res(i).scheme);
end
xlabel(ax4, 'Size (um)');
ylabel(ax4, 'Mass drop at end (%)');
title(ax4, 'Conservation behavior');
grid(ax4, 'on');

title(tl, 'Advection-only scheme comparison');
save_figure(fig, fullfile(fig_dir, 'advection_only_scheme_comparison.png'));
close(fig);
end

function log_path = save_summary_log(res, log_dir)
log_path = fullfile(log_dir, 'advection_only_scheme_comparison.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Advection-only scheme comparison\n\n');

for i = 1:numel(res)
    m = res(i).val.metrics;
    fprintf(fid, '%s\n', upper(res(i).scheme));
    fprintf(fid, '  max_cfl = %.4f\n', res(i).sim.cfl.max_cfl);
    fprintf(fid, '  min_conc = %.6e\n', m.min_conc);
    fprintf(fid, '  neg_count = %d\n', m.neg_count);
    fprintf(fid, '  mean_abs_travel_error_pct = %.3f\n', m.mean_abs_error_pct);
    fprintf(fid, '  mean_signal_width_day = %.3f\n', m.mean_signal_width_day);
    fprintf(fid, '  max_mass_drop_pct = %.3f\n\n', m.max_mass_drop_pct);
end

fprintf(fid, 'Short reading:\n');
fprintf(fid, '- lower travel error is better\n');
fprintf(fid, '- lower signal width means less artificial spreading\n');
fprintf(fid, '- low negative count and non-negative values suggest better stability\n');
fprintf(fid, '- lower mass drop means better conservation in this finite column test\n');
fclose(fid);
end

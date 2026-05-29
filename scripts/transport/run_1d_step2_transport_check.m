% run_1d_step2_transport_check
% Steps:
% 1. test the advection-only transport
% 2. compare two schemes in a simple case
% 3. save simple figures and one short summary

clear;
clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
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

law_name = 'kriest_8';
size_um = [100; 500; 1000; 3000];
pulse_amp = [1.0; 0.6; 0.3; 0.15];
size_cm = size_um .* 1e-4;
speed_cm_s = sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(speed_m_s);
cfg.t_max_s = 1.20 .* max(cfg.z_max_m ./ speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
cfg.pulse_amp = pulse_amp;
cfg.law_name = law_name;

schemes = {'upwind', 'lax_wendroff'};
scheme_labels = {'upwind', 'lax wendroff'};
plot_cols = lines(numel(schemes));
rows = struct('scheme', {}, 'size_um', {}, 'pulse_amp', {}, 'speed_m_day', {}, ...
    'expected_day', {}, 'peak_day', {}, 'error_pct', {}, 'signal_width_day', {}, ...
    'neg_count', {}, 'min_conc', {}, 'max_cfl', {});
bottom_curves = cell(numel(schemes), 1);
t_curves = cell(numel(schemes), 1);

for i = 1:numel(schemes)
    cfg.scheme = schemes{i};
    sim = solve_advection_only(cfg);
    summary = local_travel_summary(sim);

    i_rep = find(size_um == 1000, 1, 'first');
    t_curves{i} = sim.t_s ./ 86400.0;
    bottom_curves{i} = sim.bottom_signal(:, i_rep);

    neg_count = sum(sim.conc(:) < -1e-12);
    min_conc = min(sim.conc(:));
    max_cfl = sim.cfl.max_cfl;

    for j = 1:height(summary)
        row = struct();
        row.scheme = string(schemes{i});
        row.size_um = summary.size_um(j);
        row.pulse_amp = pulse_amp(j);
        row.speed_m_day = summary.speed_m_day(j);
        row.expected_day = summary.expected_day(j);
        row.peak_day = summary.peak_day(j);
        row.error_pct = summary.error_pct(j);
        row.signal_width_day = summary.signal_width_day(j);
        row.neg_count = neg_count;
        row.min_conc = min_conc;
        row.max_cfl = max_cfl;
        rows(end + 1) = row; %#ok<AGROW>
    end
end

t_max_day = 0;
for i = 1:numel(schemes)
    t_max_day = max(t_max_day, max(t_curves{i}));
end
t_ref = linspace(0, t_max_day, 2500).';
y_ref = zeros(numel(t_ref), numel(schemes));
for i = 1:numel(schemes)
    y_ref(:, i) = interp1(t_curves{i}, bottom_curves{i}, t_ref, 'linear', 0.0);
end

fig1 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax1 = axes(fig1);
hold(ax1, 'on');
for i = 1:numel(schemes)
    plot(ax1, t_ref, y_ref(:, i), 'LineWidth', 1.4, 'Color', plot_cols(i, :), ...
        'DisplayName', scheme_labels{i});
end
xlabel(ax1, 'Time (day)');
ylabel(ax1, 'Particles at 1000 m');
title(ax1, 'Advection check for 1 mm');
legend(ax1, 'Location', 'best', 'Box', 'off');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, 'step2_transport_bottom_signal_1mm.png'));
close(fig1);

T = struct2table(rows);

fig2 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax2 = axes(fig2);
hold(ax2, 'on');
for i = 1:numel(schemes)
    mask = T.scheme == schemes{i};
    plot(ax2, T.size_um(mask), abs(T.error_pct(mask)), 'LineWidth', 1.4, ...
        'Color', plot_cols(i, :), 'DisplayName', scheme_labels{i});
end
xlabel(ax2, 'Particle size (um)');
ylabel(ax2, 'Arrival time error (%)');
title(ax2, 'Advection timing error');
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, 'step2_transport_travel_error.png'));
close(fig2);

csv_path = fullfile(tab_dir, 'step2_transport_check_summary.csv');
writetable(T, csv_path);

log_path = fullfile(log_dir, 'step2_transport_check.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Step 2: advection-only transport check\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- law = %s\n', law_name);
fprintf(fid, '- sizes = 100, 500, 1000, 3000 um\n');
fprintf(fid, '- pulse amplitudes = 1.0, 0.6, 0.3, 0.15\n');
fprintf(fid, '- target CFL = 0.5\n\n');

for i = 1:numel(schemes)
    mask = T.scheme == schemes{i};
    fprintf(fid, '%s\n', upper(schemes{i}));
    fprintf(fid, '- max_cfl = %.4f\n', T.max_cfl(find(mask, 1, 'first')));
    fprintf(fid, '- min_conc = %.6e\n', T.min_conc(find(mask, 1, 'first')));
    fprintf(fid, '- neg_count = %d\n', T.neg_count(find(mask, 1, 'first')));
    fprintf(fid, '- mean abs travel-time error = %.4f %%\n', mean(abs(T.error_pct(mask))));
    fprintf(fid, '- mean signal width = %.4f day\n\n', mean(T.signal_width_day(mask)));
end

fprintf(fid, 'Reading:\n');
fprintf(fid, '- low negative count is important\n');
fprintf(fid, '- lower travel-time error is better\n');
fprintf(fid, '- lower signal width means less numerical spreading\n');
fprintf(fid, '- if one scheme gives many negative values, it is not a safe baseline yet\n');
fclose(fid);

disp('Saved step 2 figures and summary:');
disp(fullfile(fig_dir, 'step2_transport_bottom_signal_1mm.png'));
disp(fullfile(fig_dir, 'step2_transport_travel_error.png'));
disp(csv_path);
disp(log_path);

function out = local_travel_summary(sim)
t_day = sim.t_s(:) ./ 86400.0;
expected_day = sim.cfg.z_max_m ./ sim.speed_m_day;
ns = numel(sim.size_um);
peak_day = zeros(ns, 1);
signal_width_day = zeros(ns, 1);

for is = 1:ns
    y = sim.bottom_signal(:, is);
    [~, idx] = max(y);
    peak_day(is) = t_day(idx);
    signal_width_day(is) = local_width(t_day, y);
end

err_pct = 100.0 .* (peak_day - expected_day) ./ max(expected_day, realmin);
out = table(sim.size_um, sim.speed_m_day, expected_day, peak_day, err_pct, signal_width_day, ...
    'VariableNames', {'size_um', 'speed_m_day', 'expected_day', 'peak_day', 'error_pct', 'signal_width_day'});
end

function width = local_width(x, y)
w = y(:);
x = x(:);
w(~isfinite(w)) = 0;
w(w < 0) = 0;
sw = sum(w);
if sw <= 0
    width = NaN;
    return;
end
mu = sum(x .* w) ./ sw;
width = sqrt(sum(((x - mu) .^ 2) .* w) ./ sw);
end

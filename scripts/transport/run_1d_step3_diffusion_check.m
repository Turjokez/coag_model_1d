% run_1d_step3_diffusion_check
% Steps:
% 1. add diffusion on top of trusted upwind transport
% 2. compare no-diff and diff cases
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
cfg.scheme = 'upwind';
cfg.kz_m2_s = 1e-4;

sim_adv = solve_advection_only(rmfield(cfg, 'kz_m2_s'));
sim_diff = solve_advection_diffusion(cfg);

summary = local_diffusion_summary(sim_adv, sim_diff);
csv_path = fullfile(tab_dir, 'step3_diffusion_check_summary.csv');
writetable(summary, csv_path);

t_day_adv = sim_adv.t_s ./ 86400.0;
t_day_diff = sim_diff.t_s ./ 86400.0;
t_max_day = max([t_day_adv(:); t_day_diff(:)]);
t_ref = linspace(0, t_max_day, 2500).';
i_rep = find(size_um == 1000, 1, 'first');
y_adv = interp1(t_day_adv, sim_adv.bottom_signal(:, i_rep), t_ref, 'linear', 0.0);
y_diff = interp1(t_day_diff, sim_diff.bottom_signal(:, i_rep), t_ref, 'linear', 0.0);

fig1 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax1 = axes(fig1);
hold(ax1, 'on');
plot(ax1, t_ref, y_adv, 'k', 'LineWidth', 1.4, 'DisplayName', 'no diffusion');
plot(ax1, t_ref, y_diff, 'r', 'LineWidth', 1.4, 'DisplayName', 'with diffusion');
xlabel(ax1, 'Time (day)');
ylabel(ax1, 'Particles at 1000 m');
title(ax1, '1 mm signal with and without diffusion');
legend(ax1, 'Location', 'best', 'Box', 'off');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, 'step3_diffusion_bottom_signal_1mm.png'));
close(fig1);

fig2 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax2 = axes(fig2);
hold(ax2, 'on');
plot(ax2, summary.size_um, summary.base_width_day, 'k--', 'LineWidth', 1.4, 'DisplayName', 'no diffusion');
plot(ax2, summary.size_um, summary.diff_width_day, 'r', 'LineWidth', 1.4, 'DisplayName', 'with diffusion');
xlabel(ax2, 'Particle size (um)');
ylabel(ax2, 'Signal width at 1000 m (day)');
title(ax2, 'Diffusion makes the signal wider');
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, 'step3_diffusion_signal_width.png'));
close(fig2);

fig3 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax3 = axes(fig3);
hold(ax3, 'on');
adv_total = sum(sim_adv.tracked_mass, 2);
diff_total = sum(sim_diff.tracked_mass, 2);
plot(ax3, t_day_adv, tracked_error_pct(adv_total), 'k', 'LineWidth', 1.2, 'DisplayName', 'adv tracked total');
plot(ax3, t_day_diff, tracked_error_pct(diff_total), 'r', 'LineWidth', 1.2, 'DisplayName', 'diff tracked total');
plot(ax3, t_day_diff, tracked_error_pct(sim_diff.tracked_volume_total), 'b', 'LineWidth', 1.2, 'DisplayName', 'diff tracked volume');
xlabel(ax3, 'Time (day)');
ylabel(ax3, 'Volume error (%)');
title(ax3, 'Diffusion conservation');
legend(ax3, 'Location', 'best', 'Box', 'off');
ax3.LineWidth = 1.0;
ax3.FontSize = 11;
save_figure(fig3, fullfile(fig_dir, 'step3_diffusion_conservation.png'));
close(fig3);

log_path = fullfile(log_dir, 'step3_diffusion_check.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Step 3: diffusion check on trusted upwind transport\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- law = %s\n', law_name);
fprintf(fid, '- scheme = upwind\n');
fprintf(fid, '- Kz = %.3e m^2/s\n', cfg.kz_m2_s);
fprintf(fid, '- sizes = 100, 500, 1000, 3000 um\n');
fprintf(fid, '- pulse amplitudes = 1.0, 0.6, 0.3, 0.15\n');
fprintf(fid, '- diff_alpha = %.4f\n\n', sim_diff.diff_alpha);

adv_total = sum(sim_adv.tracked_mass, 2);
diff_total = sum(sim_diff.tracked_mass, 2);

fprintf(fid, 'Max tracked error:\n');
fprintf(fid, '- advection tracked total = %.6e %%\n', max(abs(tracked_error_pct(adv_total))));
fprintf(fid, '- diffusion tracked total = %.6e %%\n', max(abs(tracked_error_pct(diff_total))));
fprintf(fid, '- diffusion tracked volume = %.6e %%\n\n', max(abs(tracked_error_pct(sim_diff.tracked_volume_total))));

fprintf(fid, 'By size:\n');
for i = 1:height(summary)
    fprintf(fid, '- %g um: base width = %.4f day, diff width = %.4f day, peak shift = %.4f day\n', ...
        summary.size_um(i), summary.base_width_day(i), summary.diff_width_day(i), summary.peak_shift_day(i));
end

fprintf(fid, '\nReading:\n');
fprintf(fid, '- conservation should stay close to zero error\n');
fprintf(fid, '- with diffusion, signal width should usually be the same or larger\n');
fprintf(fid, '- if width becomes smaller, that needs another look\n');
fclose(fid);

disp('Saved step 3 figures and summary:');
disp(fullfile(fig_dir, 'step3_diffusion_bottom_signal_1mm.png'));
disp(fullfile(fig_dir, 'step3_diffusion_signal_width.png'));
disp(fullfile(fig_dir, 'step3_diffusion_conservation.png'));
disp(csv_path);
disp(log_path);

function out = local_diffusion_summary(sim_adv, sim_diff)
ns = numel(sim_adv.size_um);
t_day_adv = sim_adv.t_s ./ 86400.0;
t_day_diff = sim_diff.t_s ./ 86400.0;

base_width_day = zeros(ns, 1);
diff_width_day = zeros(ns, 1);
base_peak_day = zeros(ns, 1);
diff_peak_day = zeros(ns, 1);
peak_shift_day = zeros(ns, 1);

for is = 1:ns
    y_adv = sim_adv.bottom_signal(:, is);
    y_diff = sim_diff.bottom_signal(:, is);
    base_width_day(is) = local_width(t_day_adv, y_adv);
    diff_width_day(is) = local_width(t_day_diff, y_diff);

    [~, ia] = max(y_adv);
    [~, id] = max(y_diff);
    base_peak_day(is) = t_day_adv(ia);
    diff_peak_day(is) = t_day_diff(id);
    peak_shift_day(is) = diff_peak_day(is) - base_peak_day(is);
end

out = table(sim_adv.size_um, sim_adv.speed_m_day, base_width_day, diff_width_day, ...
    base_peak_day, diff_peak_day, peak_shift_day, ...
    'VariableNames', {'size_um', 'speed_m_day', 'base_width_day', 'diff_width_day', ...
    'base_peak_day', 'diff_peak_day', 'peak_shift_day'});
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

function err = tracked_error_pct(y)
if isvector(y)
    y = y(:);
    err = 100.0 .* (y - y(1)) ./ max(abs(y(1)), realmin);
else
    err = zeros(size(y, 1), 1);
    y0 = y(1, :);
    for i = 1:size(y, 1)
        now = y(i, :);
        now_err = 100.0 .* (now - y0) ./ max(abs(y0), realmin);
        err(i) = max(abs(now_err));
    end
end
end

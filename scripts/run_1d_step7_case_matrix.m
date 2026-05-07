% run_1d_step7_case_matrix
% Short note:
% 1. rerun step-7 depth setup with 4 process cases
% 2. isolate where large number change comes from
% 3. save simple figures + table + log

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');
tab_dir = fullfile(repo_root, 'output', 'tables');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end
if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end

% Base setup from step-7
law_name = 'kriest_8';
size_um = round(logspace(log10(200), log10(3000), 8))';
size_cm = size_um .* 1e-4;
pulse_amp = powerlaw_concentration(size_cm, 5e-3, -2.5);
base_speed_cm_s = sinking_speed_named(size_cm, law_name);
base_speed_m_s = base_speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(base_speed_m_s);
cfg.t_max_s = 1.10 .* max(cfg.z_max_m ./ base_speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = base_speed_m_s;
cfg.pulse_amp = pulse_amp;
cfg.law_name = law_name;
cfg.scheme = 'upwind';
cfg.kz_m2_s = 1e-4;
cfg.kernel_mode = 'shear_only';
cfg.epsilon_mks = 1e-6;
cfg.coag_scale = 100.0;
cfg.coag_substeps = 4;
cfg.scale_shear = 1.0;
cfg.scale_diff_sed = 0.0;
cfg.frag_substeps = 4;
cfg.c3 = 0.005;
cfg.c4 = 1.45;

grid = make_depth_grid(cfg.z_max_m, cfg.dz_m);
prof = local_profiles(grid.z_m);
sink_prof = build_sinking_speed_profile(base_speed_m_s, prof.temp_c, prof.rho_kg_m3);

cfg_depth = cfg;
cfg_depth.kz_profile_m2_s = prof.kz_m2_s;
cfg_depth.temp_profile_c = prof.temp_c;
cfg_depth.sal_profile_psu = prof.sal_psu;
cfg_depth.rho_profile_kg_m3 = prof.rho_kg_m3;
cfg_depth.speed_profile_m_s = sink_prof.speed_profile_m_s;

% 4 cases to isolate effect
case_names = {
    'transport_sink_only'
    'transport_sink_coag'
    'transport_sink_frag'
    'transport_sink_coag_frag'
};
case_labels = {
    'transport + sink'
    'transport + sink + coag'
    'transport + sink + frag'
    'transport + sink + coag + frag'
};

ncase = numel(case_names);
sim_cases = cell(ncase, 1);

for i = 1:ncase
    switch case_names{i}
        case 'transport_sink_only'
            sim_cases{i} = solve_column_1d_core(cfg_depth, true, false, false);
        case 'transport_sink_coag'
            sim_cases{i} = solve_column_1d_core(cfg_depth, true, true, false);
        case 'transport_sink_frag'
            sim_cases{i} = solve_column_1d_core(cfg_depth, true, false, true);
        case 'transport_sink_coag_frag'
            sim_cases{i} = solve_column_1d_core(cfg_depth, true, true, true);
        otherwise
            error('Unknown case: %s', case_names{i});
    end
end

t_day = sim_cases{1}.t_s(:) ./ 86400.0;
small_mask = size_um <= 500;

rows = struct('case_name', {}, 'neg_count', {}, 'max_track_err_pct', {}, ...
    'final_total_number_change_pct', {}, 'time_to_80_small_day', {});

small_curves = zeros(numel(t_day), ncase);
track_curves = zeros(numel(t_day), ncase);
number_curves = zeros(numel(t_day), ncase);

for i = 1:ncase
    sim = sim_cases{i};
    v0 = max(sim.tracked_volume_total(1), realmin);
    small_curves(:, i) = 100.0 .* sum(sim.column_volume_by_size(:, small_mask), 2) ./ v0;
    track_curves(:, i) = 100.0 .* (sim.tracked_volume_total - sim.tracked_volume_total(1)) ./ ...
        max(abs(sim.tracked_volume_total(1)), realmin);
    number_curves(:, i) = 100.0 .* sim.total_number ./ max(sim.total_number(1), realmin);

    neg_count = sum(sim.conc(:) < -1e-12);
    max_track_err_pct = max(abs(track_curves(:, i)));
    final_total_number_change_pct = 100.0 .* (sim.total_number(end) - sim.total_number(1)) ./ ...
        max(abs(sim.total_number(1)), realmin);
    time80 = local_first_cross_day(t_day, small_curves(:, i) ./ 100.0, 0.80);

    row = struct();
    row.case_name = string(case_names{i});
    row.neg_count = neg_count;
    row.max_track_err_pct = max_track_err_pct;
    row.final_total_number_change_pct = final_total_number_change_pct;
    row.time_to_80_small_day = time80;
    rows(end + 1) = row; %#ok<AGROW>
end

summary = struct2table(rows);
csv_path = fullfile(tab_dir, 'step7_case_matrix_summary.csv');
writetable(summary, csv_path);

% Figure: small-size volume
fig1 = figure('Color', 'w', 'Position', [120 120 760 540]);
ax1 = axes(fig1);
hold(ax1, 'on');
cols = lines(ncase);
for i = 1:ncase
    plot(ax1, t_day, small_curves(:, i), 'LineWidth', 1.5, ...
        'Color', cols(i, :), 'DisplayName', case_labels{i});
end
xlabel(ax1, 'Time (day)');
ylabel(ax1, 'Volume in sizes <= 500 um (%)');
title(ax1, 'Step-7 matrix: small-size volume');
legend(ax1, 'Location', 'best', 'Box', 'off');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, 'step7_case_matrix_small_size_volume.png'));
close(fig1);

% Figure: tracked-volume error
fig2 = figure('Color', 'w', 'Position', [120 120 760 540]);
ax2 = axes(fig2);
hold(ax2, 'on');
for i = 1:ncase
    plot(ax2, t_day, track_curves(:, i), 'LineWidth', 1.5, ...
        'Color', cols(i, :), 'DisplayName', case_labels{i});
end
xlabel(ax2, 'Time (day)');
ylabel(ax2, 'Tracked volume error (%)');
title(ax2, 'Step-7 matrix: conservation');
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, 'step7_case_matrix_conservation.png'));
close(fig2);

% Figure: total number (% of start)
fig3 = figure('Color', 'w', 'Position', [120 120 760 540]);
ax3 = axes(fig3);
hold(ax3, 'on');
for i = 1:ncase
    plot(ax3, t_day, number_curves(:, i), 'LineWidth', 1.5, ...
        'Color', cols(i, :), 'DisplayName', case_labels{i});
end
xlabel(ax3, 'Time (day)');
ylabel(ax3, 'Total number (% of start)');
title(ax3, 'Step-7 matrix: total-number change');
legend(ax3, 'Location', 'best', 'Box', 'off');
ax3.LineWidth = 1.0;
ax3.FontSize = 11;
save_figure(fig3, fullfile(fig_dir, 'step7_case_matrix_total_number.png'));
close(fig3);

% Log
log_path = fullfile(log_dir, 'step7_case_matrix.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Step-7 case matrix\n\n');
fprintf(fid, 'Purpose:\n');
fprintf(fid, '- isolate which process causes large number change in step-7 depth run\n');
fprintf(fid, '- keep the same depth setup and change process flags only\n\n');
fprintf(fid, 'Shared setup:\n');
fprintf(fid, '- law = %s\n', law_name);
fprintf(fid, '- depth = %.0f m, dz = %.1f m\n', cfg.z_max_m, cfg.dz_m);
fprintf(fid, '- Kz(z) active\n');
fprintf(fid, '- w(z) active via viscosity scale\n');
fprintf(fid, '- c3 = %.4f, c4 = %.2f\n\n', cfg.c3, cfg.c4);

fprintf(fid, 'Case summary:\n');
for i = 1:height(summary)
    fprintf(fid, '- %s | neg_count=%d | max_track_err=%.6e %% | final_total_number_change=%.6f %% | t80_small=%.6f day\n', ...
        char(summary.case_name(i)), summary.neg_count(i), summary.max_track_err_pct(i), ...
        summary.final_total_number_change_pct(i), summary.time_to_80_small_day(i));
end

fprintf(fid, '\nDecision rule:\n');
fprintf(fid, '- if neg_count is 0 and tracked error is small, numerics are fine\n');
fprintf(fid, '- then large total-number change comes from process setup, not solver instability\n');
fclose(fid);

disp('Saved step-7 matrix outputs:');
disp(fullfile(fig_dir, 'step7_case_matrix_small_size_volume.png'));
disp(fullfile(fig_dir, 'step7_case_matrix_conservation.png'));
disp(fullfile(fig_dir, 'step7_case_matrix_total_number.png'));
disp(csv_path);
disp(log_path);

function prof = local_profiles(z_m)
prof = struct();
prof.temp_c = 4.0 + 14.0 .* exp(-z_m ./ 150.0);
prof.sal_psu = 34.2 + 0.8 .* (1.0 - exp(-z_m ./ 300.0));
prof.rho_kg_m3 = 1024.5 + 2.2 .* (1.0 - exp(-z_m ./ 250.0));
prof.kz_m2_s = 1e-6 + (1.5e-3 - 1e-6) .* exp(-z_m ./ 120.0);
end

function out = local_first_cross_day(t_day, y, thresh)
idx = find(y >= thresh, 1, 'first');
if isempty(idx)
    out = NaN;
else
    out = t_day(idx);
end
end


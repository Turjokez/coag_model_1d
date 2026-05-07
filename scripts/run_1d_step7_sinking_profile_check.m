% run_1d_step7_sinking_profile_check
% Steps:
% 1. keep the trusted step-6 depth-profile case
% 2. add depth-dependent sinking speed
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

cfg_depth = cfg;
cfg_depth.kz_profile_m2_s = prof.kz_m2_s;
cfg_depth.temp_profile_c = prof.temp_c;
cfg_depth.sal_profile_psu = prof.sal_psu;
cfg_depth.rho_profile_kg_m3 = prof.rho_kg_m3;
sim_depth = solve_with_fragmentation(cfg_depth);

sink_prof = build_sinking_speed_profile(base_speed_m_s, prof.temp_c, prof.rho_kg_m3);
cfg_sink = cfg_depth;
cfg_sink.speed_profile_m_s = sink_prof.speed_profile_m_s;
sim_sink = solve_with_fragmentation(cfg_sink);

t_day = sim_depth.t_s ./ 86400.0;
small_mask = sim_depth.size_um <= 500;
depth_small = sum(sim_depth.column_volume_by_size(:, small_mask), 2) ./ max(sim_depth.column_volume_total, realmin);
sink_small = sum(sim_sink.column_volume_by_size(:, small_mask), 2) ./ max(sim_sink.column_volume_total, realmin);
depth_t80 = first_cross_day(t_day, depth_small, 0.80);
sink_t80 = first_cross_day(t_day, sink_small, 0.80);

summary = table(sim_depth.size_um, sim_depth.speed_m_day, ...
    sim_sink.speed_profile_m_day(1, :)', sim_sink.speed_profile_m_day(end, :)', ...
    sim_depth.column_number(end, :)', sim_sink.column_number(end, :)', ...
    'VariableNames', {'size_um', 'base_speed_m_day', 'top_speed_m_day', 'bottom_speed_m_day', ...
    'depth_final_column_number', 'sink_final_column_number'});
writetable(summary, fullfile(tab_dir, 'step7_sinking_profile_summary.csv'));

i_rep = find(size_um == 940, 1, 'first');
if isempty(i_rep)
    i_rep = find(size_um == min(size_um(size_um >= 900)), 1, 'first');
end

fig1 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax1 = axes(fig1);
plot(ax1, sink_prof.scale, grid.z_m, 'r', 'LineWidth', 1.4);
set(ax1, 'YDir', 'reverse');
xlabel(ax1, 'Speed scale');
ylabel(ax1, 'Depth (m)');
title(ax1, 'Relative sinking speed with depth');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, 'step7_sinking_scale_profile.png'));
close(fig1);

fig2 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax2 = axes(fig2);
hold(ax2, 'on');
plot(ax2, sim_depth.speed_profile_m_day(:, i_rep), grid.z_m, 'k', 'LineWidth', 1.4, 'DisplayName', 'depth profile only');
plot(ax2, sim_sink.speed_profile_m_day(:, i_rep), grid.z_m, 'r', 'LineWidth', 1.4, 'DisplayName', 'with sinking coupling');
set(ax2, 'YDir', 'reverse');
xlabel(ax2, 'Speed (m/day)');
ylabel(ax2, 'Depth (m)');
title(ax2, sprintf('Sinking speed for %.3g mm particle', size_um(i_rep) ./ 1000.0));
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, 'step7_sinking_speed_profile.png'));
close(fig2);

fig3 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax3 = axes(fig3);
hold(ax3, 'on');
plot(ax3, t_day, tracked_error_pct(sim_depth.tracked_volume_total), 'k', 'LineWidth', 1.4, 'DisplayName', 'depth profile only');
plot(ax3, t_day, tracked_error_pct(sim_sink.tracked_volume_total), 'r', 'LineWidth', 1.4, 'DisplayName', 'with sinking coupling');
xlabel(ax3, 'Time (day)');
ylabel(ax3, 'Volume error (%)');
title(ax3, 'Sinking-coupling conservation');
legend(ax3, 'Location', 'best', 'Box', 'off');
ax3.LineWidth = 1.0;
ax3.FontSize = 11;
save_figure(fig3, fullfile(fig_dir, 'step7_sinking_conservation.png'));
close(fig3);

fig4 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax4 = axes(fig4);
hold(ax4, 'on');
plot(ax4, t_day, 100.0 .* depth_small, 'k', 'LineWidth', 1.4, 'DisplayName', 'depth profile only');
plot(ax4, t_day, 100.0 .* sink_small, 'r', 'LineWidth', 1.4, 'DisplayName', 'with sinking coupling');
xlabel(ax4, 'Time (day)');
ylabel(ax4, 'Volume in sizes <= 500 um (%)');
title(ax4, 'Small-particle volume');
legend(ax4, 'Location', 'best', 'Box', 'off');
ax4.LineWidth = 1.0;
ax4.FontSize = 11;
save_figure(fig4, fullfile(fig_dir, 'step7_sinking_small_size_volume.png'));
close(fig4);

fig5 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax5 = axes(fig5);
hold(ax5, 'on');
plot(ax5, summary.size_um, summary.depth_final_column_number, 'k', 'LineWidth', 1.4, ...
    'DisplayName', 'depth profile only');
plot(ax5, summary.size_um, summary.sink_final_column_number, 'r', 'LineWidth', 1.4, ...
    'DisplayName', 'with sinking coupling');
set(ax5, 'XScale', 'log', 'YScale', 'log');
xlabel(ax5, 'Size (um)');
ylabel(ax5, 'Particles left in column');
title(ax5, 'Final size distribution');
legend(ax5, 'Location', 'best', 'Box', 'off');
ax5.LineWidth = 1.0;
ax5.FontSize = 11;
save_figure(fig5, fullfile(fig_dir, 'step7_sinking_final_column_psd.png'));
close(fig5);

fid = fopen(fullfile(log_dir, 'step7_sinking_profile.txt'), 'w');
fprintf(fid, 'Step 7: depth-dependent sinking speed\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- baseline = depth-profile case from step 6\n');
fprintf(fid, '- law_name = %s\n', law_name);
fprintf(fid, '- sinking scaling uses relative kinematic viscosity\n');
fprintf(fid, '- top speed scale = %.6f\n', sink_prof.scale(1));
fprintf(fid, '- bottom speed scale = %.6f\n', sink_prof.scale(end));
fprintf(fid, '- top nu = %.6e m^2/s\n', sink_prof.nu_m2_s(1));
fprintf(fid, '- bottom nu = %.6e m^2/s\n\n', sink_prof.nu_m2_s(end));

fprintf(fid, 'Main checks:\n');
fprintf(fid, '- max tracked-volume error, depth profile only = %.6e %%\n', max(abs(tracked_error_pct(sim_depth.tracked_volume_total))));
fprintf(fid, '- max tracked-volume error, with sinking coupling = %.6e %%\n', max(abs(tracked_error_pct(sim_sink.tracked_volume_total))));
fprintf(fid, '- time to 80%% small-size volume, depth profile only = %.6f day\n', depth_t80);
fprintf(fid, '- time to 80%% small-size volume, with sinking coupling = %.6f day\n', sink_t80);
fprintf(fid, '- final total-number change, depth profile only = %.6f %%\n', total_change_pct(sim_depth.total_number));
fprintf(fid, '- final total-number change, with sinking coupling = %.6f %%\n', total_change_pct(sim_sink.total_number));
fprintf(fid, '- final column number at 200 um, depth profile only = %.6e\n', summary.depth_final_column_number(1));
fprintf(fid, '- final column number at 200 um, with sinking coupling = %.6e\n', summary.sink_final_column_number(1));
fprintf(fid, '- final column number at 3000 um, depth profile only = %.6e\n', summary.depth_final_column_number(end));
fprintf(fid, '- final column number at 3000 um, with sinking coupling = %.6e\n\n', summary.sink_final_column_number(end));

fprintf(fid, 'Reading:\n');
fprintf(fid, '- this first coupling scales sinking speed with relative kinematic viscosity only\n');
fprintf(fid, '- colder, more viscous deep water slows sinking in the coupled run\n');
fprintf(fid, '- the effect should be stronger for the larger bins because they sink more in the first place\n');
fprintf(fid, '- use the final column PSD figure, not a raw ratio, because the baseline can be near zero in the largest bins\n');
fclose(fid);

disp('Saved step 7 figures and summary:');
disp(fullfile(fig_dir, 'step7_sinking_scale_profile.png'));
disp(fullfile(fig_dir, 'step7_sinking_speed_profile.png'));
disp(fullfile(fig_dir, 'step7_sinking_conservation.png'));
disp(fullfile(fig_dir, 'step7_sinking_small_size_volume.png'));
disp(fullfile(fig_dir, 'step7_sinking_final_column_psd.png'));
disp(fullfile(tab_dir, 'step7_sinking_profile_summary.csv'));
disp(fullfile(log_dir, 'step7_sinking_profile.txt'));

function prof = local_profiles(z_m)
prof = struct();
prof.temp_c = 4.0 + 14.0 .* exp(-z_m ./ 150.0);
prof.sal_psu = 34.2 + 0.8 .* (1.0 - exp(-z_m ./ 300.0));
prof.rho_kg_m3 = 1024.5 + 2.2 .* (1.0 - exp(-z_m ./ 250.0));
prof.kz_m2_s = 1e-6 + (1.5e-3 - 1e-6) .* exp(-z_m ./ 120.0);
end

function out = first_cross_day(t_day, y, thresh)
idx = find(y >= thresh, 1, 'first');
if isempty(idx)
    out = NaN;
else
    out = t_day(idx);
end
end

function err = tracked_error_pct(y)
y = y(:);
err = 100.0 .* (y - y(1)) ./ max(abs(y(1)), realmin);
end

function out = total_change_pct(y)
y = y(:);
out = 100.0 .* (y(end) - y(1)) ./ max(abs(y(1)), realmin);
end

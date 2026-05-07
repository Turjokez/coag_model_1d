% run_1d_step6_depth_structure_check
% Steps:
% 1. start from the trusted step-5 baseline
% 2. add simple depth profiles
% 3. use Kz(z) in transport and save the other profiles for later coupling

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
speed_cm_s = sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;

cfg = struct();
cfg.z_max_m = 1000.0;
cfg.dz_m = 5.0;
cfg.dt_s = 0.50 .* cfg.dz_m ./ max(speed_m_s);
cfg.t_max_s = 1.10 .* max(cfg.z_max_m ./ speed_m_s);
cfg.size_um = size_um;
cfg.speed_m_s = speed_m_s;
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

sim_const = solve_with_fragmentation(cfg);

grid = make_depth_grid(cfg.z_max_m, cfg.dz_m);
prof = local_profiles(grid.z_m);
cfg_depth = cfg;
cfg_depth.kz_profile_m2_s = prof.kz_m2_s;
cfg_depth.temp_profile_c = prof.temp_c;
cfg_depth.sal_profile_psu = prof.sal_psu;
cfg_depth.rho_profile_kg_m3 = prof.rho_kg_m3;
sim_depth = solve_with_fragmentation(cfg_depth);

t_day = sim_const.t_s ./ 86400.0;
small_mask = sim_const.size_um <= 500;
const_small = sum(sim_const.column_volume_by_size(:, small_mask), 2) ./ max(sim_const.column_volume_total, realmin);
depth_small = sum(sim_depth.column_volume_by_size(:, small_mask), 2) ./ max(sim_depth.column_volume_total, realmin);
const_t80 = first_cross_day(t_day, const_small, 0.80);
depth_t80 = first_cross_day(t_day, depth_small, 0.80);

summary = table(sim_const.size_um, sim_const.speed_m_day, ...
    sim_const.column_number(end, :)', sim_depth.column_number(end, :)', ...
    sim_depth.column_number(end, :)' ./ max(sim_const.column_number(end, :)', realmin), ...
    'VariableNames', {'size_um', 'speed_m_day', 'const_final_column_number', ...
    'depth_final_column_number', 'depth_to_const_ratio'});
writetable(summary, fullfile(tab_dir, 'step6_depth_structure_summary.csv'));

fig1 = figure('Color', 'w', 'Position', [100 100 1120 520]);
tl = tiledlayout(fig1, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
plot(ax1, prof.temp_c, grid.z_m, 'k', 'LineWidth', 1.4);
set(ax1, 'YDir', 'reverse');
xlabel(ax1, 'Temperature (C)');
ylabel(ax1, 'Depth (m)');
title(ax1, 'Temperature');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;

ax2 = nexttile(tl, 2);
plot(ax2, prof.sal_psu, grid.z_m, 'k', 'LineWidth', 1.4);
set(ax2, 'YDir', 'reverse');
xlabel(ax2, 'Salinity (psu)');
ylabel(ax2, 'Depth (m)');
title(ax2, 'Salinity');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;

ax3 = nexttile(tl, 3);
plot(ax3, prof.rho_kg_m3, grid.z_m, 'k', 'LineWidth', 1.4);
set(ax3, 'YDir', 'reverse');
xlabel(ax3, 'Density (kg/m^3)');
ylabel(ax3, 'Depth (m)');
title(ax3, 'Density');
ax3.LineWidth = 1.0;
ax3.FontSize = 11;

ax4 = nexttile(tl, 4);
plot(ax4, 1e3 .* prof.kz_m2_s, grid.z_m, 'r', 'LineWidth', 1.4);
set(ax4, 'YDir', 'reverse');
xlabel(ax4, 'Mixing Kz (10^-3 m^2/s)');
ylabel(ax4, 'Depth (m)');
title(ax4, 'Mixing');
ax4.LineWidth = 1.0;
ax4.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, 'step6_depth_profiles.png'));
close(fig1);

fig2 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax2 = axes(fig2);
hold(ax2, 'on');
plot(ax2, t_day, tracked_error_pct(sim_const.tracked_volume_total), 'k', 'LineWidth', 1.4, 'DisplayName', 'constant');
plot(ax2, t_day, tracked_error_pct(sim_depth.tracked_volume_total), 'r', 'LineWidth', 1.4, 'DisplayName', 'depth profile');
xlabel(ax2, 'Time (day)');
ylabel(ax2, 'Volume error (%)');
title(ax2, 'Depth-profile conservation');
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, 'step6_depth_conservation.png'));
close(fig2);

fig3 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax3 = axes(fig3);
hold(ax3, 'on');
plot(ax3, t_day, 100.0 .* const_small, 'k', 'LineWidth', 1.4, 'DisplayName', 'constant');
plot(ax3, t_day, 100.0 .* depth_small, 'r', 'LineWidth', 1.4, 'DisplayName', 'depth profile');
xlabel(ax3, 'Time (day)');
ylabel(ax3, 'Volume in sizes <= 500 um (%)');
title(ax3, 'Small-particle volume');
legend(ax3, 'Location', 'best', 'Box', 'off');
ax3.LineWidth = 1.0;
ax3.FontSize = 11;
save_figure(fig3, fullfile(fig_dir, 'step6_depth_small_size_volume.png'));
close(fig3);

fig4 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax4 = axes(fig4);
hold(ax4, 'on');
plot(ax4, summary.size_um, max(sim_const.column_number(end, :)', realmin), 'k', ...
    'LineWidth', 1.4, 'DisplayName', 'constant water column');
plot(ax4, summary.size_um, max(sim_depth.column_number(end, :)', realmin), 'r', ...
    'LineWidth', 1.4, 'DisplayName', 'depth-dependent column');
set(ax4, 'XScale', 'log', 'YScale', 'log');
xlabel(ax4, 'Particle size (um)');
ylabel(ax4, 'Particles left in column');
title(ax4, 'Final size distribution');
legend(ax4, 'Location', 'best', 'Box', 'off');
ax4.LineWidth = 1.0;
ax4.FontSize = 11;
save_figure(fig4, fullfile(fig_dir, 'step6_depth_column_psd.png'));
close(fig4);

fid = fopen(fullfile(log_dir, 'step6_depth_structure.txt'), 'w');
fprintf(fid, 'Step 6: depth-dependent water-column structure\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- baseline solver = upwind + diffusion + coag + frag\n');
fprintf(fid, '- law_name = %s\n', law_name);
fprintf(fid, '- kernel_mode = shear_only\n');
fprintf(fid, '- constant Kz baseline = %.3e m^2/s\n', cfg.kz_m2_s);
fprintf(fid, '- depth-profile Kz range = %.3e to %.3e m^2/s\n', min(prof.kz_m2_s), max(prof.kz_m2_s));
fprintf(fid, '- temp range = %.2f to %.2f C\n', min(prof.temp_c), max(prof.temp_c));
fprintf(fid, '- sal range = %.2f to %.2f psu\n', min(prof.sal_psu), max(prof.sal_psu));
fprintf(fid, '- rho range = %.2f to %.2f kg/m^3\n\n', min(prof.rho_kg_m3), max(prof.rho_kg_m3));

fprintf(fid, 'Main checks:\n');
fprintf(fid, '- max tracked-volume error, constant = %.6e %%\n', max(abs(tracked_error_pct(sim_const.tracked_volume_total))));
fprintf(fid, '- max tracked-volume error, depth profile = %.6e %%\n', max(abs(tracked_error_pct(sim_depth.tracked_volume_total))));
fprintf(fid, '- time to 80%% small-size volume, constant = %.6f day\n', const_t80);
fprintf(fid, '- time to 80%% small-size volume, depth profile = %.6f day\n', depth_t80);
fprintf(fid, '- final total-number change, constant = %.6f %%\n', total_change_pct(sim_const.total_number));
fprintf(fid, '- final total-number change, depth profile = %.6f %%\n\n', total_change_pct(sim_depth.total_number));

fprintf(fid, 'Reading:\n');
fprintf(fid, '- in this first step, only Kz(z) is coupled into the transport\n');
fprintf(fid, '- T(z), S(z), and rho(z) are now saved as depth-dependent model fields\n');
fprintf(fid, '- later work can couple those profiles into settling or kernel physics if needed\n');
fclose(fid);

disp('Saved step 6 figures and summary:');
disp(fullfile(fig_dir, 'step6_depth_profiles.png'));
disp(fullfile(fig_dir, 'step6_depth_conservation.png'));
disp(fullfile(fig_dir, 'step6_depth_small_size_volume.png'));
disp(fullfile(fig_dir, 'step6_depth_column_psd.png'));
disp(fullfile(tab_dir, 'step6_depth_structure_summary.csv'));
disp(fullfile(log_dir, 'step6_depth_structure.txt'));

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

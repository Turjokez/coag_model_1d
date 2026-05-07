% run_scheme_compare_steps3to5
% Short note:
% 1. compare upwind and lax_wendroff for step 3 to 5
% 2. keep one simple figure only
% 3. use trust checks, not too many metrics

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

rows = {};

rows(end + 1, :) = local_run_diffusion('upwind'); %#ok<SAGROW>
rows(end + 1, :) = local_run_diffusion('lax_wendroff'); %#ok<SAGROW>
rows(end + 1, :) = local_run_coagulation('upwind'); %#ok<SAGROW>
rows(end + 1, :) = local_run_coagulation('lax_wendroff'); %#ok<SAGROW>
rows(end + 1, :) = local_run_fragmentation('upwind'); %#ok<SAGROW>
rows(end + 1, :) = local_run_fragmentation('lax_wendroff'); %#ok<SAGROW>

summary = cell2table(rows, 'VariableNames', ...
    {'stage', 'scheme', 'neg_count', 'min_conc', 'max_track_err_pct'});

csv_path = fullfile(tab_dir, 'scheme_compare_steps3to5.csv');
writetable(summary, csv_path);

stages = {'diffusion', 'coagulation', 'fragmentation'};
up_idx = strcmp(summary.scheme, 'upwind');
lw_idx = strcmp(summary.scheme, 'lax_wendroff');

up_neg = local_pick(summary.neg_count, summary.stage, summary.scheme, stages, 'upwind');
lw_neg = local_pick(summary.neg_count, summary.stage, summary.scheme, stages, 'lax_wendroff');
up_min = local_pick(summary.min_conc, summary.stage, summary.scheme, stages, 'upwind');
lw_min = local_pick(summary.min_conc, summary.stage, summary.scheme, stages, 'lax_wendroff');
up_err = local_pick(summary.max_track_err_pct, summary.stage, summary.scheme, stages, 'upwind');
lw_err = local_pick(summary.max_track_err_pct, summary.stage, summary.scheme, stages, 'lax_wendroff');

fig = figure('Color', 'w', 'Position', [90 90 1200 420]);
tl = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
local_group_bar(ax1, up_neg, lw_neg, stages);
set(ax1, 'YScale', 'log');
ylabel(ax1, 'Negative value count');
title(ax1, 'Negativity check');
legend(ax1, {'upwind', 'lax\_wendroff'}, 'Location', 'northwest', 'Box', 'off');

ax2 = nexttile(tl, 2);
local_group_bar(ax2, up_min, lw_min, stages);
ylabel(ax2, 'Minimum concentration');
title(ax2, 'Most negative value');

ax3 = nexttile(tl, 3);
local_group_bar(ax3, up_err, lw_err, stages);
ylabel(ax3, 'Max tracked-volume error (%)');
title(ax3, 'Conservation check');

title(tl, 'Scheme compare for step 3 to 5');
save_figure(fig, fullfile(fig_dir, 'scheme_compare_steps3to5.png'));
close(fig);

log_path = fullfile(log_dir, 'scheme_compare_steps3to5.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Scheme compare for diffusion, coagulation, and fragmentation\n\n');
fprintf(fid, 'Reading note:\n');
fprintf(fid, '- negative value count should stay at zero if possible\n');
fprintf(fid, '- minimum concentration shows how far below zero the run went\n');
fprintf(fid, '- tracked-volume error should stay close to zero\n\n');
for i = 1:height(summary)
    fprintf(fid, '%s | %s | neg_count = %d | min_conc = %.6e | max_track_err_pct = %.6e\n', ...
        summary.stage{i}, summary.scheme{i}, summary.neg_count(i), summary.min_conc(i), summary.max_track_err_pct(i));
end
fclose(fid);

disp('Saved scheme compare outputs:');
disp(fullfile(fig_dir, 'scheme_compare_steps3to5.png'));
disp(csv_path);
disp(log_path);

function row = local_run_diffusion(scheme)
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
cfg.scheme = scheme;
cfg.kz_m2_s = 1e-4;

sim = solve_advection_diffusion(cfg);
row = {char("diffusion"), char(string(scheme)), ...
    sum(sim.conc(:) < -1e-12), min(sim.conc(:)), local_max_err(sim.tracked_volume_total)};
end

function row = local_run_coagulation(scheme)
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
cfg.scheme = scheme;
cfg.kz_m2_s = 1e-4;
cfg.kernel_mode = 'shear_only';
cfg.epsilon_mks = 1e-6;
cfg.coag_scale = 100.0;
cfg.coag_substeps = 4;
cfg.scale_shear = 1.0;
cfg.scale_diff_sed = 0.0;

sim = solve_with_coagulation(cfg);
row = {char("coagulation"), char(string(scheme)), ...
    sum(sim.conc(:) < -1e-12), min(sim.conc(:)), local_max_err(sim.tracked_volume_total)};
end

function row = local_run_fragmentation(scheme)
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
cfg.scheme = scheme;
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

sim = solve_with_fragmentation(cfg);
row = {char("fragmentation"), char(string(scheme)), ...
    sum(sim.conc(:) < -1e-12), min(sim.conc(:)), local_max_err(sim.tracked_volume_total)};
end

function out = local_max_err(y)
y = y(:);
out = max(abs(100.0 .* (y - y(1)) ./ max(abs(y(1)), realmin)));
end

function vals = local_pick(x, stage_col, scheme_col, stage_names, want_scheme)
vals = zeros(numel(stage_names), 1);
for i = 1:numel(stage_names)
    idx = strcmp(stage_col, stage_names{i}) & strcmp(scheme_col, want_scheme);
    vals(i) = x(idx);
end
end

function local_group_bar(ax, y1, y2, labels)
x = 1:numel(labels);
bar(ax, x - 0.16, y1, 0.32, 'FaceColor', [0.00 0.45 0.74], 'EdgeColor', 'none');
hold(ax, 'on');
bar(ax, x + 0.16, y2, 0.32, 'FaceColor', [0.85 0.33 0.10], 'EdgeColor', 'none');
set(ax, 'XTick', x, 'XTickLabel', labels);
ax.LineWidth = 1.0;
ax.FontSize = 11;
grid(ax, 'on');
box(ax, 'on');
end

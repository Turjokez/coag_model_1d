% run_apr29_step4_kernel_mode_compare
% Short note:
% 1. rerun the same Step 4 case with three kernels
% 2. compare budget, final PSD, number, and size-volume split
% 3. keep this as a 1-D test only

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
size_um = round(logspace(log10(200), log10(3000), 8))';
size_cm = size_um .* 1e-4;
pulse_amp = local_powerlaw_concentration(size_cm, 5e-3, -2.5);
speed_cm_s = local_sinking_speed_named(size_cm, law_name);
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

sim_base = local_solve_advection_diffusion(cfg);

mode_names = {'shear_only', 'diff_sed_only', 'shear_plus_diff_sed'};
mode_labels = {'shear only', 'diff sed only', 'shear + diff sed'};
mode_cols = [0.85 0.10 0.10; 0.10 0.35 0.85; 0.10 0.55 0.20];
sim_modes = cell(numel(mode_names), 1);

for i = 1:numel(mode_names)
    cfg_now = cfg;
    cfg_now.kernel_mode = mode_names{i};
    cfg_now.epsilon_mks = 1e-6;
    cfg_now.coag_scale = 100.0;
    cfg_now.coag_substeps = 4;
    cfg_now.scale_shear = 1.0;
    cfg_now.scale_diff_sed = 1.0;
    sim_modes{i} = local_solve_with_coagulation(cfg_now);
end

t_day = sim_base.t_s ./ 86400.0;
init_volume = sim_base.tracked_volume_total(1);
init_number = sim_base.total_number(1);
small_mask = size_um < 500;
large_mask = size_um >= 500;

rows = struct('kernel_mode', {}, 'max_tracked_volume_error_pct', {}, ...
    'neg_count', {}, 'min_conc', {}, 'beta_max_m3_s', {}, ...
    'final_total_number_change_pct', {}, 'final_small_volume_pct_init', {}, ...
    'final_large_volume_pct_init', {}, 'final_tracked_volume_pct_init', {});

row = struct();
row.kernel_mode = string('no_coag');
row.max_tracked_volume_error_pct = max(abs(percent_change(sim_base.tracked_volume_total)));
row.neg_count = sum(sim_base.conc(:) < -1e-12);
row.min_conc = min(sim_base.conc(:));
row.beta_max_m3_s = NaN;
row.final_total_number_change_pct = percent_change(sim_base.total_number(end), sim_base.total_number(1));
row.final_small_volume_pct_init = 100.0 .* sum(sim_base.column_volume_by_size(end, small_mask), 2) ./ init_volume;
row.final_large_volume_pct_init = 100.0 .* sum(sim_base.column_volume_by_size(end, large_mask), 2) ./ init_volume;
row.final_tracked_volume_pct_init = 100.0 .* sim_base.tracked_volume_total(end) ./ init_volume;
rows(end + 1) = row;

for i = 1:numel(mode_names)
    sim = sim_modes{i};
    row = struct();
    row.kernel_mode = string(mode_names{i});
    row.max_tracked_volume_error_pct = max(abs(percent_change(sim.tracked_volume_total)));
    row.neg_count = sum(sim.conc(:) < -1e-12);
    row.min_conc = min(sim.conc(:));
    row.beta_max_m3_s = max(sim.beta_m3_s(:));
    row.final_total_number_change_pct = percent_change(sim.total_number(end), sim.total_number(1));
    row.final_small_volume_pct_init = 100.0 .* sum(sim.column_volume_by_size(end, small_mask), 2) ./ init_volume;
    row.final_large_volume_pct_init = 100.0 .* sum(sim.column_volume_by_size(end, large_mask), 2) ./ init_volume;
    row.final_tracked_volume_pct_init = 100.0 .* sim.tracked_volume_total(end) ./ init_volume;
    rows(end + 1) = row; %#ok<SAGROW>
end

T = struct2table(rows);
csv_path = fullfile(tab_dir, 'apr29_step4_kernel_mode_compare.csv');
writetable(T, csv_path);

make_budget_figure(fig_dir, t_day, init_volume, sim_modes, mode_labels);
make_psd_figure(fig_dir, size_um, sim_base, sim_modes, mode_labels, mode_cols);
make_number_figure(fig_dir, t_day, init_number, sim_base, sim_modes, mode_labels, mode_cols);
make_size_volume_figure(fig_dir, t_day, init_volume, small_mask, large_mask, sim_base, sim_modes, mode_labels, mode_cols);

log_path = fullfile(log_dir, 'apr29_step4_kernel_mode_compare.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Step 4 kernel-mode comparison\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- same setup as the Step 4 coagulation check\n');
fprintf(fid, '- law = %s\n', law_name);
fprintf(fid, '- scheme = upwind\n');
fprintf(fid, '- Kz = %.1e m2/s\n', cfg.kz_m2_s);
fprintf(fid, '- epsilon = %.1e m2/s3\n', 1e-6);
fprintf(fid, '- coag_scale = %.1f\n', 100.0);
fprintf(fid, '- modes = shear_only, diff_sed_only, shear_plus_diff_sed\n\n');

fprintf(fid, 'Main results:\n');
for i = 1:height(T)
    fprintf(fid, '%s\n', char(T.kernel_mode(i)));
    fprintf(fid, '- max tracked-volume error = %.6e %%\n', T.max_tracked_volume_error_pct(i));
    fprintf(fid, '- neg_count = %d\n', T.neg_count(i));
    fprintf(fid, '- beta max = %.6e m3/s\n', T.beta_max_m3_s(i));
    fprintf(fid, '- final total-number change = %.6f %%\n', T.final_total_number_change_pct(i));
    fprintf(fid, '- final small volume = %.6f %% of initial total volume\n', T.final_small_volume_pct_init(i));
    fprintf(fid, '- final large volume = %.6f %% of initial total volume\n\n', T.final_large_volume_pct_init(i));
end

fprintf(fid, 'Reading:\n');
fprintf(fid, '- this is a matched comparison, not a tuned comparison\n');
fprintf(fid, '- all three kernel modes use the same transport, grid, initial PSD, Kz, and coag_scale\n');
fprintf(fid, '- if tracked volume stays near 100%%, the kernel mode is not breaking conservation\n');
fprintf(fid, '- the main physical comparison is the final PSD and number loss relative to the no-coag baseline\n');
fprintf(fid, '- the mixed kernel should be read after checking the separate shear-only and diff-sed-only cases\n');
fclose(fid);

disp('Saved Step 4 kernel-mode comparison:');
disp(fullfile(fig_dir, 'apr29_step4_kernel_budget.png'));
disp(fullfile(fig_dir, 'apr29_step4_kernel_final_psd.png'));
disp(fullfile(fig_dir, 'apr29_step4_kernel_total_number.png'));
disp(fullfile(fig_dir, 'apr29_step4_kernel_size_volume.png'));
disp(csv_path);
disp(log_path);

function make_budget_figure(fig_dir, t_day, init_volume, sim_modes, mode_labels)
fig = figure('Color', 'w', 'Position', [120 120 1100 430]);
tl = tiledlayout(fig, 1, numel(sim_modes), 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(sim_modes)
    sim = sim_modes{i};
    ax = nexttile(tl, i);
    hold(ax, 'on');
    plot(ax, t_day, 100.0 .* sim.column_volume_total ./ init_volume, ...
        'k', 'LineWidth', 1.3, 'DisplayName', 'volume in column');
    plot(ax, t_day, 100.0 .* sim.export_volume_total ./ init_volume, ...
        'r', 'LineWidth', 1.3, 'DisplayName', 'volume out bottom');
    plot(ax, t_day, 100.0 .* sim.tracked_volume_total ./ init_volume, ...
        'b', 'LineWidth', 1.3, 'DisplayName', 'tracked total');
    xlabel(ax, 'Time (day)');
    ylabel(ax, 'Volume (% initial)');
    title(ax, mode_labels{i});
    ylim(ax, [0, 105]);
    if i == 1
        legend(ax, 'Location', 'best', 'Box', 'off');
    end
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
end
title(tl, 'Budget and conservation');
local_save_figure(fig, fullfile(fig_dir, 'apr29_step4_kernel_budget.png'));
close(fig);
end

function make_psd_figure(fig_dir, size_um, sim_base, sim_modes, mode_labels, mode_cols)
all_num = [sim_base.column_number(1, :), sim_base.column_number(end, :)];
for i = 1:numel(sim_modes)
    all_num = [all_num, sim_modes{i}.column_number(end, :)]; %#ok<AGROW>
end
floor_val = 1e-12 .* max([all_num(:); realmin]);

fig = figure('Color', 'w', 'Position', [120 120 760 560]);
ax = axes(fig);
hold(ax, 'on');
plot(ax, size_um, max(sim_base.column_number(1, :), floor_val), ...
    'k--', 'LineWidth', 1.2, 'DisplayName', 'initial');
plot(ax, size_um, max(sim_base.column_number(end, :), floor_val), ...
    'k', 'LineWidth', 1.4, 'DisplayName', 'no coag');
for i = 1:numel(sim_modes)
    plot(ax, size_um, max(sim_modes{i}.column_number(end, :), floor_val), ...
        'LineWidth', 1.4, 'Color', mode_cols(i, :), 'DisplayName', mode_labels{i});
end
plot(ax, size_um, floor_val .* ones(size(size_um)), ':', ...
    'Color', [0.45 0.45 0.45], 'LineWidth', 1.0, 'DisplayName', 'plot floor');
set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, 'Particle size (um)');
ylabel(ax, 'Particles left in column');
title(ax, 'Final size distribution');
legend(ax, 'Location', 'southwest', 'Box', 'off');
grid(ax, 'on');
ax.LineWidth = 1.0;
ax.FontSize = 11;
local_save_figure(fig, fullfile(fig_dir, 'apr29_step4_kernel_final_psd.png'));
close(fig);
end

function make_number_figure(fig_dir, t_day, init_number, sim_base, sim_modes, mode_labels, mode_cols)
fig = figure('Color', 'w', 'Position', [120 120 720 520]);
ax = axes(fig);
hold(ax, 'on');
plot(ax, t_day, 100.0 .* sim_base.total_number ./ init_number, ...
    'k', 'LineWidth', 1.4, 'DisplayName', 'no coag');
for i = 1:numel(sim_modes)
    plot(ax, t_day, 100.0 .* sim_modes{i}.total_number ./ init_number, ...
        'LineWidth', 1.4, 'Color', mode_cols(i, :), 'DisplayName', mode_labels{i});
end
xlabel(ax, 'Time (day)');
ylabel(ax, 'Total number (% initial)');
title(ax, 'Total number left in column');
legend(ax, 'Location', 'best', 'Box', 'off');
ax.LineWidth = 1.0;
ax.FontSize = 11;
local_save_figure(fig, fullfile(fig_dir, 'apr29_step4_kernel_total_number.png'));
close(fig);
end

function make_size_volume_figure(fig_dir, t_day, init_volume, small_mask, large_mask, sim_base, sim_modes, mode_labels, mode_cols)
fig = figure('Color', 'w', 'Position', [120 120 1080 460]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
hold(ax1, 'on');
plot(ax1, t_day, 100.0 .* sum(sim_base.column_volume_by_size(:, small_mask), 2) ./ init_volume, ...
    'k', 'LineWidth', 1.4, 'DisplayName', 'no coag');
for i = 1:numel(sim_modes)
    plot(ax1, t_day, 100.0 .* sum(sim_modes{i}.column_volume_by_size(:, small_mask), 2) ./ init_volume, ...
        'LineWidth', 1.4, 'Color', mode_cols(i, :), 'DisplayName', mode_labels{i});
end
xlabel(ax1, 'Time (day)');
ylabel(ax1, 'Volume (% initial total)');
title(ax1, 'Small sizes < 500 um');
legend(ax1, 'Location', 'best', 'Box', 'off');

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
plot(ax2, t_day, 100.0 .* sum(sim_base.column_volume_by_size(:, large_mask), 2) ./ init_volume, ...
    'k', 'LineWidth', 1.4, 'DisplayName', 'no coag');
for i = 1:numel(sim_modes)
    plot(ax2, t_day, 100.0 .* sum(sim_modes{i}.column_volume_by_size(:, large_mask), 2) ./ init_volume, ...
        'LineWidth', 1.4, 'Color', mode_cols(i, :), 'DisplayName', mode_labels{i});
end
xlabel(ax2, 'Time (day)');
ylabel(ax2, 'Volume (% initial total)');
title(ax2, 'Large sizes >= 500 um');

title(tl, 'Column volume by size group');
local_save_figure(fig, fullfile(fig_dir, 'apr29_step4_kernel_size_volume.png'));
close(fig);
end

function out = percent_change(y, y0)
if nargin < 2
    y0 = y(1);
end
out = 100.0 .* (y - y0) ./ max(abs(y0), realmin);
end

function sim = local_solve_advection_diffusion(cfg)
sim = local_solve_column_1d(cfg, false);
end

function sim = local_solve_with_coagulation(cfg)
sim = local_solve_column_1d(cfg, true);
end

function sim = local_solve_column_1d(cfg, do_coag)
z_m = (0:cfg.dz_m:cfg.z_max_m)';
nz = numel(z_m);
t_s = (0:cfg.dt_s:cfg.t_max_s)';
nt = numel(t_s);

size_um = cfg.size_um(:);
ns = numel(size_um);
speed_m_s = cfg.speed_m_s(:);
pulse_amp = cfg.pulse_amp(:)';

conc = zeros(nz, ns, nt);
top_mask = z_m <= 50.0;
conc(top_mask, :, 1) = repmat(pulse_amp, sum(top_mask), 1);

d_m = size_um .* 1e-6;
vol_part_m3 = (pi/6) .* (d_m .^ 3);

column_number = zeros(nt, ns);
column_volume_by_size = zeros(nt, ns);
export_volume_total = zeros(nt, 1);
export_number = zeros(ns, 1);

beta_m3_s = zeros(ns, ns);
if do_coag
    beta_m3_s = local_build_beta_matrix(size_um .* 1e-4, cfg);
end

for is = 1:ns
    col_num = sum(conc(:, is, 1)) .* cfg.dz_m;
    column_number(1, is) = col_num;
    column_volume_by_size(1, is) = col_num .* vol_part_m3(is);
end
column_volume_total = zeros(nt, 1);
tracked_volume_total = zeros(nt, 1);
total_number = zeros(nt, 1);
column_volume_total(1) = sum(column_volume_by_size(1, :));
tracked_volume_total(1) = column_volume_total(1);
total_number(1) = sum(column_number(1, :));

kz = 0.0;
if isfield(cfg, 'kz_m2_s') && ~isempty(cfg.kz_m2_s)
    kz = cfg.kz_m2_s;
end
alpha_d = kz .* cfg.dt_s ./ max(cfg.dz_m .* cfg.dz_m, realmin);

n_sub = 1;
if isfield(cfg, 'coag_substeps') && ~isempty(cfg.coag_substeps)
    n_sub = max(1, round(cfg.coag_substeps));
end
dt_sub = cfg.dt_s ./ n_sub;

for it = 2:nt
    c_now = conc(:, :, it - 1);
    c_new = zeros(nz, ns);

    for is = 1:ns
        c_col = local_upwind_flux_step(c_now(:, is), speed_m_s(is), cfg.dt_s, cfg.dz_m, 0.0);

        if alpha_d > 0
            c_tmp = c_col;
            c_tmp(2:end-1) = c_col(2:end-1) + alpha_d .* ...
                (c_col(3:end) - 2 .* c_col(2:end-1) + c_col(1:end-2));
            c_tmp(1) = c_tmp(2);
            c_tmp(end) = c_tmp(end-1);
            c_col = c_tmp;
        end

        c_col(c_col < 0) = 0;
        c_new(:, is) = c_col;

        bottom_flux = max(speed_m_s(is), 0.0) .* c_now(end, is);
        export_number(is) = export_number(is) + bottom_flux .* cfg.dt_s;
    end

    if do_coag
        for iz = 1:nz
            n_vec = c_new(iz, :)';
            for sub = 1:n_sub
                dn = zeros(ns, 1);
                for i = 1:ns
                    for j = i:ns
                        if i == j
                            coll = beta_m3_s(i, j) .* n_vec(i) .* n_vec(j) .* dt_sub;
                            coll = min(coll, 0.25 .* n_vec(i));
                            k = min(ns, i + 1);
                            dn(i) = dn(i) - 2.0 .* coll;
                            dn(k) = dn(k) + coll;
                        else
                            coll = beta_m3_s(i, j) .* n_vec(i) .* n_vec(j) .* dt_sub;
                            coll = min(coll, 0.25 .* min(n_vec(i), n_vec(j)));
                            k = min(ns, max(i, j) + 1);
                            dn(i) = dn(i) - coll;
                            dn(j) = dn(j) - coll;
                            dn(k) = dn(k) + coll;
                        end
                    end
                end
                n_vec = n_vec + dn;
                n_vec(n_vec < 0) = 0;
            end
            c_new(iz, :) = n_vec';
        end
    end

    conc(:, :, it) = c_new;

    for is = 1:ns
        col_num = sum(c_new(:, is)) .* cfg.dz_m;
        column_number(it, is) = col_num;
        column_volume_by_size(it, is) = col_num .* vol_part_m3(is);
    end

    export_volume_total(it) = sum(export_number .* vol_part_m3);
    column_volume_total(it) = sum(column_volume_by_size(it, :));
    tracked_volume_total(it) = column_volume_total(it) + export_volume_total(it);
    total_number(it) = sum(column_number(it, :));
end

sim = struct();
sim.t_s = t_s;
sim.z_m = z_m;
sim.size_um = size_um;
sim.conc = conc;
sim.column_number = column_number;
sim.total_number = total_number;
sim.column_volume_by_size = column_volume_by_size;
sim.column_volume_total = column_volume_total;
sim.export_volume_total = export_volume_total;
sim.tracked_volume_total = tracked_volume_total;
sim.beta_m3_s = beta_m3_s;
end

function beta_m3_s = local_build_beta_matrix(size_cm, cfg)
[D1, D2] = ndgrid(size_cm, size_cm);

mode_name = "shear_only";
if isfield(cfg, 'kernel_mode') && ~isempty(cfg.kernel_mode)
    mode_name = lower(string(cfg.kernel_mode));
end

[beta_ds_cm3_s, ~, ~] = local_beta_diff_sed_from_law(D1, D2, cfg.law_name);
beta_ds_m3_s = beta_ds_cm3_s .* 1e-6;

eps_mks = 1e-6;
if isfield(cfg, 'epsilon_mks') && ~isempty(cfg.epsilon_mks)
    eps_mks = cfg.epsilon_mks;
end
rg_m = 0.5 .* (D1 + D2) .* 1e-2;
beta_shear_m3_s = sqrt(max(eps_mks, 0)) .* (rg_m .^ 3);

switch mode_name
    case "shear_only"
        beta_m3_s = beta_shear_m3_s;
    case "diff_sed_only"
        beta_m3_s = beta_ds_m3_s;
    otherwise
        beta_m3_s = beta_shear_m3_s + beta_ds_m3_s;
end

if isfield(cfg, 'scale_shear') && ~isempty(cfg.scale_shear)
    if mode_name ~= "diff_sed_only"
        beta_m3_s = beta_m3_s .* cfg.scale_shear;
    end
end
if isfield(cfg, 'scale_diff_sed') && ~isempty(cfg.scale_diff_sed)
    if mode_name ~= "shear_only"
        beta_m3_s = beta_m3_s .* cfg.scale_diff_sed;
    end
end
if isfield(cfg, 'coag_scale') && ~isempty(cfg.coag_scale)
    beta_m3_s = beta_m3_s .* cfg.coag_scale;
end

beta_m3_s(~isfinite(beta_m3_s)) = 0;
beta_m3_s(beta_m3_s < 0) = 0;
end

function [beta, w1, w2] = local_beta_diff_sed_from_law(d1_cm, d2_cm, law_name)
w1 = local_sinking_speed_named(d1_cm, law_name);
w2 = local_sinking_speed_named(d2_cm, law_name);
beta = (pi/4) .* (d1_cm + d2_cm) .* (d1_cm + d2_cm) .* abs(w1 - w2);
end

function c = local_powerlaw_concentration(size_cm, c0, slope)
if nargin < 2 || isempty(c0)
    c0 = 1.0;
end
if nargin < 3 || isempty(slope)
    slope = -2.5;
end
ref = min(size_cm(size_cm > 0));
if isempty(ref)
    ref = 1.0;
end
c = c0 .* (size_cm ./ ref) .^ slope;
c(~isfinite(c)) = 0;
c(c < 0) = 0;
end

function w_cm_s = local_sinking_speed_named(d_cm, law_name)
cfg = SimulationConfig();
law = lower(string(law_name));
switch law
    case "current"
        w_cm_s = local_current_law(d_cm, cfg);
    case "kriest_8"
        w_cm_s = (66 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "kriest_9"
        w_cm_s = (132 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "siegel_2025"
        d_mm = d_cm .* 10.0;
        w_cm_s = (20.2 .* (d_mm .^ 0.67) .* 100) ./ cfg.day_to_sec;
    otherwise
        error('Unknown law: %s', char(law_name));
end
w_cm_s(~isfinite(w_cm_s)) = 0;
w_cm_s(w_cm_s < 0) = 0;
end

function w_cm_s = local_current_law(d_cm, cfg)
r_v = 0.5 .* d_cm;
setcon = KernelLibrary.currentSetcon(cfg);
r_i = KernelLibrary.conservativeToFractalRadius(r_v, cfg);
w_cm_s = setcon .* (r_v .^ 3) ./ max(r_i, realmin);
end

function c_next = local_upwind_flux_step(c_prev, w_m_s, dt_s, dz_m, c_in)
if nargin < 5
    c_in = 0.0;
end
c_prev = c_prev(:);
n = numel(c_prev);
if isscalar(w_m_s)
    w_col = w_m_s .* ones(n, 1);
else
    w_col = w_m_s(:);
end
flux = zeros(n + 1, 1);
flux(1) = max(w_col(1), 0.0) .* c_in;
for j = 1:(n - 1)
    w_face = 0.5 .* (w_col(j) + w_col(j + 1));
    if w_face >= 0
        c_up = c_prev(j);
    else
        c_up = c_prev(j + 1);
    end
    flux(j + 1) = w_face .* c_up;
end
flux(n + 1) = max(w_col(end), 0.0) .* c_prev(end);
c_next = c_prev - (dt_s / dz_m) .* (flux(2:end) - flux(1:end-1));
c_next(c_next < 0) = 0;
end

function local_save_figure(fig_handle, fig_path)
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

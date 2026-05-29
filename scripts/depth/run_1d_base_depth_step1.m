% run_1d_base_depth_step1
% Short note:
% 1. run one clean 1-D depth-dependent base case
% 2. check pulse sinking and conservation
% 3. save simple figures, table, and pass/fail log

clear;
clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
tab_dir = fullfile(repo_root, 'output', 'tables');
log_dir = fullfile(repo_root, 'output', 'logs');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(tab_dir, 'dir'), mkdir(tab_dir); end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end

% -------- setup --------
law_name = 'kriest_8';
z_max_m = 1000.0;
dz_m = 5.0;
z_m = (0:dz_m:z_max_m).';
nz = numel(z_m);

size_um = [100; 300; 1000; 3000];
ns = numel(size_um);
size_cm = size_um .* 1e-4;
base_speed_cm_s = local_sinking_speed_named(size_cm, law_name);
base_speed_m_s = base_speed_cm_s .* 0.01;

% viscosity-based depth scaling (simple first coupling)
prof = local_profiles(z_m);
nu_rel = prof.nu_m2_s ./ max(prof.nu_m2_s(1), realmin);
speed_scale = 1.0 ./ max(nu_rel, realmin);

% speed profile per size
speed_profile_m_s = zeros(nz, ns);
for is = 1:ns
    speed_profile_m_s(:, is) = base_speed_m_s(is) .* speed_scale;
end

dt_adv_s = 0.45 .* dz_m ./ max(speed_profile_m_s(:));
dt_diff_s = 0.20 .* (dz_m .* dz_m) ./ max(prof.kz_m2_s);
dt_s = min(dt_adv_s, dt_diff_s);
t_max_day = 70.0;
t_s = (0:dt_s:(t_max_day .* 86400.0)).';
t_day = t_s ./ 86400.0;
nt = numel(t_s);

% top pulse
top_layer_m = 50.0;
pulse_amp = [1.0; 0.7; 0.4; 0.2];
conc = zeros(nz, ns, nt);
top_mask = z_m <= top_layer_m;
conc(top_mask, :, 1) = repmat(reshape(pulse_amp, 1, []), sum(top_mask), 1);

% for tracked-volume conservation
d_m = size_um .* 1e-6;
vol_part_m3 = (pi/6) .* (d_m .^ 3);
export_number = zeros(ns, 1);
column_volume_total = zeros(nt, 1);
export_volume_total = zeros(nt, 1);
tracked_volume_total = zeros(nt, 1);
neg_count = 0;
for is = 1:ns
    ncol = sum(conc(:, is, 1)) .* dz_m;
    column_volume_total(1) = column_volume_total(1) + ncol .* vol_part_m3(is);
end
tracked_volume_total(1) = column_volume_total(1);

for it = 2:nt
    c_prev = conc(:, :, it - 1);
    c_next = zeros(nz, ns);

    for is = 1:ns
        c_col = c_prev(:, is);
        w_col = speed_profile_m_s(:, is);

        % advection flux (upwind)
        adv_flux = local_adv_flux_upwind(c_col, w_col, 0.0);

        % diffusion flux (variable Kz)
        diff_flux = local_diff_flux(c_col, prof.kz_m2_s, dz_m);

        % conservative update
        c_new = c_col - (dt_s ./ dz_m) .* (adv_flux(2:end) - adv_flux(1:end-1)) ...
                      - (dt_s ./ dz_m) .* (diff_flux(2:end) - diff_flux(1:end-1));

        % bottom export bookkeeping from old state
        export_number(is) = export_number(is) + max(w_col(end), 0.0) .* c_col(end) .* dt_s;

        neg_big = c_new < -1e-12;
        neg_count = neg_count + sum(neg_big);
        c_new(c_new < 0) = 0;
        c_next(:, is) = c_new;
    end

    conc(:, :, it) = c_next;

    col_v = 0.0;
    for is = 1:ns
        ncol = sum(c_next(:, is)) .* dz_m;
        col_v = col_v + ncol .* vol_part_m3(is);
    end
    column_volume_total(it) = col_v;
    export_volume_total(it) = sum(export_number .* vol_part_m3);
    tracked_volume_total(it) = column_volume_total(it) + export_volume_total(it);
end

% -------- checks --------
rep_size_um = 1000;
[~, i_rep] = min(abs(size_um - rep_size_um));

z_center = zeros(nt, 1);
for it = 1:nt
    y = conc(:, i_rep, it);
    sw = sum(y);
    if sw > 0
        z_center(it) = sum(z_m .* y) ./ sw;
    else
        z_center(it) = NaN;
    end
end

valid = isfinite(z_center);
dzc = diff(z_center(valid));
is_monotone_downward = all(dzc >= -1e-8);

v0 = tracked_volume_total(1);
track_err_pct = 100.0 .* (tracked_volume_total - v0) ./ max(abs(v0), realmin);
max_abs_track_err_pct = max(abs(track_err_pct));

pass_neg = (neg_count == 0);
pass_track = (max_abs_track_err_pct < 1e-6);
pass_center = is_monotone_downward;
all_pass = pass_neg && pass_track && pass_center;

% -------- outputs --------
% Fig 1: profile movement for 1 mm
profile_times_day = [0, 10, 25, 45, 60];
profile_idx = local_nearest_time_idx(t_day, profile_times_day);
fig1 = figure('Color', 'w', 'Position', [120 120 640 540]);
ax1 = axes(fig1);
hold(ax1, 'on');
cols = lines(numel(profile_idx));
for i = 1:numel(profile_idx)
    it = profile_idx(i);
    plot(ax1, conc(:, i_rep, it), z_m, 'LineWidth', 1.4, ...
        'Color', cols(i, :), 'DisplayName', sprintf('%.0f day', t_day(it)));
end
set(ax1, 'YDir', 'reverse');
xlabel(ax1, 'Concentration');
ylabel(ax1, 'Depth (m)');
title(ax1, sprintf('Pulse move | %.0f um | %s', size_um(i_rep), strrep(law_name, '_', ' ')));
legend(ax1, 'Location', 'southeast', 'Box', 'off');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
local_save_figure(fig1, fullfile(fig_dir, 'base_depth_step1_pulse_profiles.png'));
close(fig1);

% Fig 2: depth-size snapshots
snap_times_day = [0, 10, 30, 60];
snap_idx = local_nearest_time_idx(t_day, snap_times_day);
all_vals = conc(:, :, snap_idx);
cmax = max(all_vals(:));
fig2 = figure('Color', 'w', 'Position', [120 120 900 700]);
tl = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(snap_idx)
    ax = nexttile(tl, i);
    imagesc(ax, 1:ns, z_m, conc(:, :, snap_idx(i)));
    axis(ax, 'xy');
    set(ax, 'YDir', 'reverse');
    caxis(ax, [0, cmax]);
    xticks(ax, 1:ns);
    xticklabels(ax, string(size_um));
    xlabel(ax, 'Size (um)');
    ylabel(ax, 'Depth (m)');
    title(ax, sprintf('t = %.0f day', t_day(snap_idx(i))));
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
end
colormap(fig2, parula);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Conc';
title(tl, 'Depth-size snapshots');
local_save_figure(fig2, fullfile(fig_dir, 'base_depth_step1_depth_size_snapshots.png'));
close(fig2);

% Fig 3: tracked volume error
fig3 = figure('Color', 'w', 'Position', [120 120 640 500]);
ax3 = axes(fig3);
plot(ax3, t_day, track_err_pct, 'k', 'LineWidth', 1.4);
xlabel(ax3, 'Time (day)');
ylabel(ax3, 'Tracked volume error (%)');
title(ax3, 'Conservation check');
ax3.LineWidth = 1.0;
ax3.FontSize = 11;
local_save_figure(fig3, fullfile(fig_dir, 'base_depth_step1_conservation.png'));
close(fig3);

% table
T = table(size_um, base_speed_m_s .* 86400.0, ...
    speed_profile_m_s(1, :)' .* 86400.0, ...
    speed_profile_m_s(end, :)' .* 86400.0, ...
    'VariableNames', {'size_um', 'base_speed_m_day', 'top_speed_m_day', 'bottom_speed_m_day'});
writetable(T, fullfile(tab_dir, 'base_depth_step1_speed_summary.csv'));

% log
log_path = fullfile(log_dir, 'base_depth_step1_check.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Base depth step 1 check\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- law = %s\n', law_name);
fprintf(fid, '- depth = %.0f m, dz = %.1f m\n', z_max_m, dz_m);
fprintf(fid, '- advection = upwind conservative flux\n');
fprintf(fid, '- diffusion = variable Kz(z) flux form\n');
fprintf(fid, '- sinking = depth-scaled with relative viscosity\n');
fprintf(fid, '- max CFL = %.4f\n\n', max(speed_profile_m_s(:)) .* dt_s ./ dz_m);
fprintf(fid, '- dt_adv_s = %.3f s\n', dt_adv_s);
fprintf(fid, '- dt_diff_s = %.3f s\n\n', dt_diff_s);

fprintf(fid, 'Pass checks:\n');
fprintf(fid, '- no negative concentration: %d (neg_count=%d)\n', pass_neg, neg_count);
fprintf(fid, '- tracked volume error < 1e-6 %%: %d (max=%.6e %%)\n', pass_track, max_abs_track_err_pct);
fprintf(fid, '- pulse center moves down monotonically: %d\n\n', pass_center);

if all_pass
    fprintf(fid, 'Final status: PASS\n');
else
    fprintf(fid, 'Final status: NEEDS WORK\n');
end
fclose(fid);

disp('Saved base depth step 1 outputs:');
disp(fullfile(fig_dir, 'base_depth_step1_pulse_profiles.png'));
disp(fullfile(fig_dir, 'base_depth_step1_depth_size_snapshots.png'));
disp(fullfile(fig_dir, 'base_depth_step1_conservation.png'));
disp(fullfile(tab_dir, 'base_depth_step1_speed_summary.csv'));
disp(log_path);

% ---------- local functions ----------
function idx = local_nearest_time_idx(t_day, target_day)
idx = zeros(size(target_day));
for i = 1:numel(target_day)
    [~, idx(i)] = min(abs(t_day - target_day(i)));
end
end

function prof = local_profiles(z_m)
prof = struct();
prof.temp_c = 4.0 + 14.0 .* exp(-z_m ./ 150.0);
prof.sal_psu = 34.2 + 0.8 .* (1.0 - exp(-z_m ./ 300.0));
prof.rho_kg_m3 = 1024.5 + 2.2 .* (1.0 - exp(-z_m ./ 250.0));
prof.kz_m2_s = 1e-6 + (1.5e-3 - 1e-6) .* exp(-z_m ./ 120.0);

% simple viscosity fit with T and S
Tk = prof.temp_c + 273.15;
mu_w = 2.414e-5 .* 10 .^ (247.8 ./ (Tk - 140.0));
mu = mu_w .* (1.0 + 0.0015 .* (prof.sal_psu - 35.0));
prof.nu_m2_s = mu ./ prof.rho_kg_m3;
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
a0 = cfg.d0 / 2.0;
amfrac_temp = (4.0/3.0*pi)^(-1.0/cfg.fr_dim) * a0^(1.0 - 3.0/cfg.fr_dim);
amfrac = amfrac_temp * sqrt(0.6);
bmfrac = 1.0 / cfg.fr_dim;
del_rho = (4.5*2.48) * cfg.kvisc * cfg.rho_fl / cfg.g * (cfg.d0/2.0)^(-0.83);
setcon = (2.0/9.0) * del_rho / cfg.rho_fl * cfg.g / cfg.kvisc;
r_v = 0.5 .* d_cm;
av_vol = (4.0/3.0) .* pi .* (r_v .^ 3);
r_i = amfrac .* (av_vol .^ bmfrac);
w_cm_s = setcon .* (r_v .^ 3) ./ max(r_i, realmin);
end

function flux = local_adv_flux_upwind(c_prev, w_m_s, c_in)
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
end

function flux = local_diff_flux(c_prev, kz_m2_s, dz_m)
c_prev = c_prev(:);
kz_m2_s = kz_m2_s(:);
n = numel(c_prev);
flux = zeros(n + 1, 1);
flux(1) = 0.0;
for j = 1:(n - 1)
    k_face = 0.5 .* (kz_m2_s(j) + kz_m2_s(j + 1));
    grad = (c_prev(j + 1) - c_prev(j)) ./ dz_m;
    flux(j + 1) = -k_face .* grad;
end
flux(n + 1) = 0.0;
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

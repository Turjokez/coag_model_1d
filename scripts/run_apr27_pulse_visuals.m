% run_apr27_pulse_visuals
% Short note:
% 1. make the pulse plots from the meeting
% 2. keep this as a visual check, not a new model result
% 3. save figures only inside the 1-D testing folder

clear;
clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
tab_dir = fullfile(repo_root, 'output', 'tables');
log_dir = fullfile(repo_root, 'output', 'logs');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(tab_dir, 'dir')
    mkdir(tab_dir);
end
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end

law_name = 'kriest_8';
law_label = strrep(law_name, '_', ' ');
z_max_m = 1000.0;
dz_m = 5.0;
z_m = (0:dz_m:z_max_m).';
nz = numel(z_m);

size_um = [100; 300; 1000; 3000];
size_cm = size_um .* 1e-4;
speed_cm_s = local_sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;
speed_m_day = speed_m_s .* 86400.0;

dt_s = 0.45 .* dz_m ./ max(speed_m_s);
t_max_day = 70.0;
t_s = (0:dt_s:(t_max_day .* 86400.0)).';
t_day = t_s ./ 86400.0;
nt = numel(t_s);

top_layer_m = 50.0;
pulse_amp = [1.0; 0.7; 0.4; 0.2];
conc = zeros(nz, numel(size_um), nt);
top_mask = z_m <= top_layer_m;
conc(top_mask, :, 1) = repmat(reshape(pulse_amp, 1, []), sum(top_mask), 1);

for it = 2:nt
    for is = 1:numel(size_um)
        speed_col = speed_m_s(is) .* ones(nz, 1);
        conc(:, is, it) = local_upwind_flux_step(conc(:, is, it - 1), speed_col, dt_s, dz_m, 0.0);
    end
end

profile_size_um = 1000;
[~, i_profile] = min(abs(size_um - profile_size_um));
profile_times_day = [0, 10, 25, 45, 60];
profile_idx = nearest_time_idx(t_day, profile_times_day);

fig1 = figure('Color', 'w', 'Position', [120 120 650 560]);
ax1 = axes(fig1);
hold(ax1, 'on');
cols = lines(numel(profile_idx));
for i = 1:numel(profile_idx)
    it = profile_idx(i);
    plot(ax1, conc(:, i_profile, it), z_m, 'LineWidth', 1.4, ...
        'Color', cols(i, :), 'DisplayName', sprintf('%.0f day', t_day(it)));
end
set(ax1, 'YDir', 'reverse');
xlabel(ax1, 'Concentration');
ylabel(ax1, 'Depth (m)');
title(ax1, sprintf('Pulse movement | %.0f um | %s', size_um(i_profile), law_label));
legend(ax1, 'Location', 'southeast', 'Box', 'off');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;
local_save_figure(fig1, fullfile(fig_dir, 'apr27_pulse_profiles.png'));
close(fig1);

snap_times_day = [0, 10, 30, 60];
snap_idx = nearest_time_idx(t_day, snap_times_day);
all_vals = conc(:, :, snap_idx);
cmax = max(all_vals(:));

fig2 = figure('Color', 'w', 'Position', [120 120 900 700]);
tl = tiledlayout(fig2, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(snap_idx)
    ax = nexttile(tl, i);
    imagesc(ax, 1:numel(size_um), z_m, conc(:, :, snap_idx(i)));
    axis(ax, 'xy');
    set(ax, 'YDir', 'reverse');
    caxis(ax, [0, cmax]);
    xticks(ax, 1:numel(size_um));
    xticklabels(ax, string(size_um));
    xlabel(ax, 'Particle size (um)');
    ylabel(ax, 'Depth (m)');
    title(ax, sprintf('t = %.0f day', t_day(snap_idx(i))));
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
end
colormap(fig2, parula);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Concentration';
title(tl, sprintf('Depth-size pulse snapshots | %s', law_label));
local_save_figure(fig2, fullfile(fig_dir, 'apr27_depth_size_snapshots.png'));
close(fig2);

T = table(size_um, speed_m_day, 1000.0 ./ speed_m_day, ...
    'VariableNames', {'size_um', 'speed_m_day', 'time_to_1000_day'});
writetable(T, fullfile(tab_dir, 'apr27_pulse_speed_summary.csv'));

log_path = fullfile(log_dir, 'apr27_pulse_visuals.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'April 27 pulse visual check\n\n');
fprintf(fid, 'Purpose:\n');
fprintf(fid, '- show a pulse moving down the 1-D column\n');
fprintf(fid, '- show depth-size snapshots at a few times\n');
fprintf(fid, '- keep this separate from the main model code\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- law = %s\n', law_name);
fprintf(fid, '- scheme = upwind flux form\n');
fprintf(fid, '- no diffusion, no coagulation, no fragmentation\n');
fprintf(fid, '- depth = %.0f m, dz = %.1f m\n', z_max_m, dz_m);
fprintf(fid, '- top pulse depth = %.0f m\n', top_layer_m);
fprintf(fid, '- max CFL = %.4f\n\n', max(speed_m_s) .* dt_s ./ dz_m);
fprintf(fid, 'Files:\n');
fprintf(fid, '- output/figures/apr27_pulse_profiles.png\n');
fprintf(fid, '- output/figures/apr27_depth_size_snapshots.png\n');
fprintf(fid, '- output/tables/apr27_pulse_speed_summary.csv\n');
fclose(fid);

disp('Saved pulse visual check:');
disp(fullfile(fig_dir, 'apr27_pulse_profiles.png'));
disp(fullfile(fig_dir, 'apr27_depth_size_snapshots.png'));
disp(fullfile(tab_dir, 'apr27_pulse_speed_summary.csv'));
disp(log_path);

function idx = nearest_time_idx(t_day, target_day)
idx = zeros(size(target_day));
for i = 1:numel(target_day)
    [~, idx(i)] = min(abs(t_day - target_day(i)));
end
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

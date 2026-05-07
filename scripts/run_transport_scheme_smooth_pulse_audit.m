% run_transport_scheme_smooth_pulse_audit
% Short note:
% 1. use one smooth pulse away from the boundaries
% 2. compare upwind and Lax-Wendroff against the exact shift
% 3. keep the test simple so transport quality is easy to read

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
size_um = [100; 500; 1000; 3000];
size_cm = size_um .* 1e-4;
speed_cm_s = sinking_speed_named(size_cm, law_name);
speed_m_s = speed_cm_s .* 0.01;
speed_m_day = speed_m_s .* 86400.0;

z_max_m = 1000.0;
dz_m = 5.0;
z_m = (0:dz_m:z_max_m).';
nz = numel(z_m);

dt_s = 0.50 .* dz_m ./ max(speed_m_s);
t_end_day = 5.0;
t_end_s = t_end_day .* 86400.0;
nt = floor(t_end_s ./ dt_s);
t_end_s = nt .* dt_s;
t_end_day = t_end_s ./ 86400.0;

center0_m = 250.0;
sigma_m = 20.0;
u0 = exp(-0.5 .* ((z_m - center0_m) ./ sigma_m) .^ 2);

schemes = {'upwind', 'lax_wendroff'};
scheme_labels = {'upwind', 'lax wendroff'};
plot_cols = lines(numel(schemes));

rows = struct('scheme', {}, 'size_um', {}, 'speed_m_day', {}, ...
    'l2_error', {}, 'linf_error', {}, 'min_conc', {}, ...
    'neg_count', {}, 'tracked_mass_error_pct', {});

rep_size_um = 1000;
i_rep = find(size_um == rep_size_um, 1, 'first');
if isempty(i_rep)
    i_rep = 1;
end

rep_profiles = struct();

for i = 1:numel(schemes)
    step_fun = pick_step_function(schemes{i});

    for is = 1:numel(size_um)
        alpha = speed_m_s(is) .* dt_s ./ dz_m;
        c_now = u0;
        tracked0 = dz_m .* sum(c_now);

        for it = 1:nt
            c_now = step_fun(c_now, alpha, 0.0);
        end

        shift_m = speed_m_s(is) .* t_end_s;
        exact_now = exp(-0.5 .* ((z_m - (center0_m + shift_m)) ./ sigma_m) .^ 2);

        err = c_now - exact_now;
        l2_error = sqrt(mean(err .^ 2));
        linf_error = max(abs(err));
        min_conc = min(c_now);
        neg_count = sum(c_now < -1e-12);
        tracked_now = dz_m .* sum(c_now);
        tracked_mass_error_pct = 100.0 .* (tracked_now - tracked0) ./ max(tracked0, realmin);

        row = struct();
        row.scheme = string(schemes{i});
        row.size_um = size_um(is);
        row.speed_m_day = speed_m_day(is);
        row.l2_error = l2_error;
        row.linf_error = linf_error;
        row.min_conc = min_conc;
        row.neg_count = neg_count;
        row.tracked_mass_error_pct = tracked_mass_error_pct;
        rows(end + 1) = row; %#ok<AGROW>

        if is == i_rep
            rep_profiles(i).scheme = string(schemes{i}); %#ok<SAGROW>
            rep_profiles(i).num = c_now; %#ok<SAGROW>
            rep_profiles(i).exact = exact_now; %#ok<SAGROW>
        end
    end
end

T = struct2table(rows);
csv_path = fullfile(tab_dir, 'transport_scheme_smooth_pulse_audit.csv');
writetable(T, csv_path);

fig1 = figure('Color', 'w', 'Position', [120 120 1200 520]);
tl = tiledlayout(fig1, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tl, 1);
hold(ax1, 'on');
plot(ax1, z_m, rep_profiles(1).exact, 'k--', 'LineWidth', 1.4, 'DisplayName', 'exact');
for i = 1:numel(schemes)
    plot(ax1, z_m, rep_profiles(i).num, 'LineWidth', 1.4, ...
        'Color', plot_cols(i, :), 'DisplayName', scheme_labels{i});
end
xlabel(ax1, 'Depth (m)');
ylabel(ax1, 'Concentration');
title(ax1, sprintf('Smooth pulse after %.2f day | %.1f mm', t_end_day, size_um(i_rep) ./ 1000.0));
legend(ax1, 'Location', 'best', 'Box', 'off');
ax1.LineWidth = 1.0;
ax1.FontSize = 11;

ax2 = nexttile(tl, 2);
hold(ax2, 'on');
for i = 1:numel(schemes)
    mask = T.scheme == schemes{i};
    plot(ax2, T.size_um(mask), T.l2_error(mask), 'LineWidth', 1.4, ...
        'Color', plot_cols(i, :), 'DisplayName', scheme_labels{i});
end
xlabel(ax2, 'Particle size (um)');
ylabel(ax2, 'L2 error');
title(ax2, 'Error against exact shifted pulse');
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0;
ax2.FontSize = 11;

save_figure(fig1, fullfile(fig_dir, 'transport_scheme_smooth_pulse_audit.png'));
close(fig1);

log_path = fullfile(log_dir, 'transport_scheme_smooth_pulse_audit.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'Smooth-pulse transport audit\n\n');
fprintf(fid, 'Setup:\n');
fprintf(fid, '- law = %s\n', law_name);
fprintf(fid, '- sizes = 100, 500, 1000, 3000 um\n');
fprintf(fid, '- smooth Gaussian pulse\n');
fprintf(fid, '- pulse center = %.1f m\n', center0_m);
fprintf(fid, '- pulse sigma = %.1f m\n', sigma_m);
fprintf(fid, '- no inflow\n');
fprintf(fid, '- no boundary hit during the run\n');
fprintf(fid, '- target CFL = 0.5\n');
fprintf(fid, '- run time = %.4f day\n\n', t_end_day);

for i = 1:numel(schemes)
    mask = T.scheme == schemes{i};
    fprintf(fid, '%s\n', upper(schemes{i}));
    fprintf(fid, '- mean L2 error = %.6e\n', mean(T.l2_error(mask)));
    fprintf(fid, '- max Linf error = %.6e\n', max(T.linf_error(mask)));
    fprintf(fid, '- min concentration = %.6e\n', min(T.min_conc(mask)));
    fprintf(fid, '- total neg_count = %d\n', sum(T.neg_count(mask)));
    fprintf(fid, '- max tracked-mass error = %.6e %%\n\n', max(abs(T.tracked_mass_error_pct(mask))));
end

fprintf(fid, 'Reading:\n');
fprintf(fid, '- this test removes the sharp top pulse and the boundary exit\n');
fprintf(fid, '- it is only checking the transport stencil itself\n');
fprintf(fid, '- if Lax-Wendroff is implemented well, it should usually beat upwind on this smooth case\n');
fclose(fid);

disp('Saved smooth-pulse audit files:');
disp(fullfile(fig_dir, 'transport_scheme_smooth_pulse_audit.png'));
disp(csv_path);
disp(log_path);

function step_fun = pick_step_function(scheme)
switch lower(string(scheme))
    case "upwind"
        step_fun = @upwind_step;
    case "lax_wendroff"
        step_fun = @lax_wendroff_step;
    otherwise
        error('run_transport_scheme_smooth_pulse_audit:scheme', ...
            'Unknown scheme: %s', scheme);
end
end

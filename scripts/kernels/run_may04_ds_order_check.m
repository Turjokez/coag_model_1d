% run_may04_ds_order_check
% Short note:
% 1. verify DS ordering for b = 0.5, 1.0, 2.0
% 2. use same prefactor so only exponent changes
% 3. save one table, one log, one figure

clear;
clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
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

b_list = [0.5, 1.0, 2.0];
d1_um = 1.0;
d2_um = [10, 100, 1000, 10000];
d1_cm = d1_um * 1e-4;
d2_cm = d2_um * 1e-4;

rows = struct('d1_um', {}, 'd2_um', {}, ...
    'beta_b0p5', {}, 'beta_b1p0', {}, 'beta_b2p0', {}, 'is_order_ok', {});

for j = 1:numel(d2_um)
    beta_vals = zeros(size(b_list));
    for i = 1:numel(b_list)
        b = b_list(i);
        w1 = d1_um .^ b;
        w2 = d2_um(j) .^ b;
        beta_vals(i) = local_beta_diff_sed(d1_cm, d2_cm(j), w1, w2);
    end

    row = struct();
    row.d1_um = d1_um;
    row.d2_um = d2_um(j);
    row.beta_b0p5 = beta_vals(1);
    row.beta_b1p0 = beta_vals(2);
    row.beta_b2p0 = beta_vals(3);
    row.is_order_ok = (beta_vals(3) > beta_vals(2)) && (beta_vals(2) > beta_vals(1));
    rows(end + 1) = row; %#ok<AGROW>
end

T = struct2table(rows);
csv_path = fullfile(tab_dir, 'may04_ds_order_check.csv');
writetable(T, csv_path);

fig = figure('Color', 'w', 'Position', [120 120 760 520]);
ax = axes(fig);
hold(ax, 'on');
plot_cols = lines(numel(b_list));
for i = 1:numel(b_list)
    b = b_list(i);
    w1 = d1_um .^ b;
    w2 = d2_um .^ b;
    beta = local_beta_diff_sed(d1_cm, d2_cm, w1, w2);
    plot(ax, d2_um, beta, 'LineWidth', 1.5, 'Color', plot_cols(i, :), ...
        'DisplayName', sprintf('b = %.1f', b));
end
set(ax, 'XScale', 'log', 'YScale', 'log');
xlabel(ax, 'Partner size d2 (um)');
ylabel(ax, 'beta ds (arb.)');
title(ax, 'DS ordering check with same prefactor');
legend(ax, 'Location', 'northwest', 'Box', 'off');
grid(ax, 'on');
save_figure(fig, fullfile(fig_dir, 'may04_ds_order_check.png'));
close(fig);

ok_all = all(T.is_order_ok);
log_path = fullfile(log_dir, 'may04_ds_order_check.txt');
fid = fopen(log_path, 'w');
fprintf(fid, 'May 04 DS order check\n\n');
fprintf(fid, 'Rule to pass:\n');
fprintf(fid, '- beta(b=2.0) > beta(b=1.0) > beta(b=0.5)\n');
fprintf(fid, '- checked at d1 = %.1f um and d2 = [10, 100, 1000, 10000] um\n\n', d1_um);

for j = 1:height(T)
    fprintf(fid, 'd2=%g um | b0.5=%.6e b1.0=%.6e b2.0=%.6e | order_ok=%d\n', ...
        T.d2_um(j), T.beta_b0p5(j), T.beta_b1p0(j), T.beta_b2p0(j), T.is_order_ok(j));
end

fprintf(fid, '\nOverall pass: %d\n', ok_all);
fclose(fid);

disp('Saved DS order check:');
disp(fullfile(fig_dir, 'may04_ds_order_check.png'));
disp(csv_path);
disp(log_path);

function beta = local_beta_diff_sed(d1_cm, d2_cm, w1, w2)
beta = (pi / 4.0) .* (d1_cm + d2_cm) .* (d1_cm + d2_cm) .* abs(w1 - w2);
end

function save_figure(fig_handle, fig_path)
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

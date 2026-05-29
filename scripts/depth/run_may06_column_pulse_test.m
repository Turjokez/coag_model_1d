% run_may06_column_pulse_test
% Basic check for the 1-D ColumnSimulation class.
%
% Tests:
%   1. Pulse starts at top, sinks down over 60 days.
%   2. No negative concentrations.
%   3. Total biovolume conserved (pulse has not reached bottom in 60 days).
%   4. Pulse center moves deeper over time.

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end

% ---- config ----
cfg = SimulationConfig();
cfg.n_sections     = 5;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.enable_coag    = false;
cfg.enable_disagg  = false;
cfg.enable_sinking = true;
cfg.enable_linear  = false;
cfg.r_to_rg        = 1.6;
cfg.t_init         = 0;
cfg.t_final        = 60;
cfg.delta_t        = 1;

% ---- depth grid ----
H   = 1000;   % m
n_z = 20;
cgrid = ColumnGrid(H, n_z);

% ---- ocean profile ----
prof = DepthProfile.typical(cgrid.z_centers);

% ---- run ----
sim    = ColumnSimulation(cfg, cgrid, prof);
result = sim.run();

Y_hist = result.concentrations;   % n_t x n_z x n_sec
t      = result.time;
n_t    = length(t);
n_z    = cgrid.n_z;

% ---- checks ----
fprintf('\n=== Pulse test checks ===\n');

% 1. negatives
neg_cnt = sum(Y_hist(:) < -1e-30);
fprintf('neg_count:        %d   (expected 0)\n', neg_cnt);

% 2. biovolume change (should decrease as particles leave bottom)
av_vol = sim.size_grid.av_vol(:)';
biovol = zeros(n_t, 1);
for it = 1:n_t
    Y_t = squeeze(Y_hist(it, :, :));
    biovol(it) = sum(Y_t .* repmat(av_vol, n_z, 1), 'all');
end
biovol_change_pct = (biovol(end) - biovol(1)) / max(biovol(1), 1e-60) * 100;
fprintf('biovolume change: %.3f %%  (expected <= 0; zero ok if pulse has not reached bottom)\n', biovol_change_pct);

% 3. pulse center depth for middle bin
bin_track = round(cfg.n_sections / 2);
for it = [1, floor(n_t/2), n_t]
    depth_prof = squeeze(Y_hist(it, :, bin_track));
    [~, iz]    = max(depth_prof);
    fprintf('bin %d center at t=%.0f d: depth = %.0f m\n', ...
        bin_track, t(it), cgrid.z_centers(iz));
end

% 4. CFL
fprintf('CFL:              %.4f\n', result.cfl);

% ---- figure: depth profiles at several times ----
t_plot  = [0, 20, 40, 60];
colors  = lines(length(t_plot));
bin_fig = round(cfg.n_sections / 2);

figure;
hold on;
for ip = 1:length(t_plot)
    [~, it] = min(abs(t - t_plot(ip)));
    prof_plot = squeeze(Y_hist(it, :, bin_fig));
    plot(prof_plot / max(prof_plot + 1e-60), cgrid.z_centers, ...
        'LineWidth', 1.8, 'Color', colors(ip, :), ...
        'DisplayName', sprintf('t = %d d', t(it)));
end
set(gca, 'YDir', 'reverse');
xlabel('Normalized concentration');
ylabel('Depth [m]');
title(sprintf('Pulse sinking: bin %d of %d  (kriest\\_8)', bin_fig, cfg.n_sections));
legend('Location', 'southeast');

fname = fullfile(fig_dir, 'may06_column_pulse_test.png');
saveas(gcf, fname);
fprintf('\nFigure saved: %s\n', fname);

% run_may11_scaling_test
% Test run time vs number of size sections.

clear; close all; clc;

% Put classes on path.
addpath('src');
repo_root = pwd;
if ~exist('SimulationConfig', 'class')
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(fullfile(repo_root, 'src')));
end

% Section counts to test.
n_sec_list = [5, 10, 15, 20, 25];
time_sec   = zeros(size(n_sec_list));

for i = 1:numel(n_sec_list)
    n_sec = n_sec_list(i);

    cfg = SimulationConfig( ...
        'n_sections', n_sec, ...
        't_final', 60, ...
        'delta_t', 1, ...
        'sinking_law', 'kriest_8', ...
        'ds_kernel_mode', 'sinking_law', ...
        'enable_coag', true, ...
        'enable_disagg', false, ...
        'proc_substeps', 20);

    col_grid = ColumnGrid(1000, 20);

    % Use the depth coordinate used by this repo.
    if isprop(col_grid, 'z_mid')
        profile = DepthProfile.typical(col_grid.z_mid);
    else
        profile = DepthProfile.typical(col_grid.z_centers);
    end

    sim = ColumnSimulation(cfg, col_grid, profile);

    t0 = tic;
    sim.run();
    time_sec(i) = toc(t0);
end

% Print table.
fprintf('%10s  %12s\n', 'n_sections', 'time_sec');
for i = 1:numel(n_sec_list)
    fprintf('%10d  %12.2f\n', n_sec_list(i), time_sec(i));
end

% Save simple scaling figure.
fig_dir = fullfile(repo_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

figure;
plot(n_sec_list, time_sec, 'o-');
xlabel('n sections');
ylabel('run time (s)');
title('scaling test — 60 days');
saveas(gcf, fullfile(fig_dir, 'may11_scaling_test.png'));


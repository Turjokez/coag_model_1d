% run_1d_zoo_production
% 1-D column test with surface production and zooplankton grazing.

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z_plot   = col_grid.z_centers;

cfg_base = SimulationConfig( ...
    'n_sections', 20, ...
    't_final', 180, ...
    'delta_t', 1, ...
    'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_coag', true, ...
    'enable_sinking', true, ...
    'enable_disagg', false, ...
    'proc_substeps', 20, ...
    'enable_surface_pp', true, ...
    'surface_pp_bin', 1, ...
    'surface_pp_rate', 3e-8);

% case 1: production only
cfg1 = cfg_base.copy();
cfg1.enable_zoo = false;
sim1 = ColumnSimulation(cfg1, col_grid, profile);
out1 = sim1.run();
Y1   = squeeze(out1.concentrations);
Y1_0 = squeeze(Y1(1, :, :));
Y1_f = squeeze(Y1(end, :, :));

% case 2: production + grazing (Zc=100, Zf=50)
cfg2 = cfg_base.copy();
cfg2.enable_zoo = true;
cfg2.zoo_Zc = 100;
cfg2.zoo_Zf = 50;
cfg2.zoo_p  = 0.3;
cfg2.zoo_ic = 1;
cfg2.zoo_c  = 1e-4;
cfg2.zoo_s  = 1e-4;
sim2 = ColumnSimulation(cfg2, col_grid, profile);
out2 = sim2.run();
Y2   = squeeze(out2.concentrations);
Y2_0 = squeeze(Y2(1, :, :));
Y2_f = squeeze(Y2(end, :, :));

% case 3: production + grazing (Zc=200, Zf=100)
cfg3 = cfg_base.copy();
cfg3.enable_zoo = true;
cfg3.zoo_Zc = 200;
cfg3.zoo_Zf = 100;
cfg3.zoo_p  = 0.3;
cfg3.zoo_ic = 1;
cfg3.zoo_c  = 1e-4;
cfg3.zoo_s  = 1e-4;
sim3 = ColumnSimulation(cfg3, col_grid, profile);
out3 = sim3.run();
Y3   = squeeze(out3.concentrations);
Y3_0 = squeeze(Y3(1, :, :));
Y3_f = squeeze(Y3(end, :, :));

% totals
bv1_0 = sum(Y1_0, 'all'); bv1_f = sum(Y1_f, 'all');
bv2_0 = sum(Y2_0, 'all'); bv2_f = sum(Y2_f, 'all');
bv3_0 = sum(Y3_0, 'all'); bv3_f = sum(Y3_f, 'all');

fprintf('Case 1 no grazing   : t=0 %.4e, t=180 %.4e\n', bv1_0, bv1_f);
fprintf('Case 2 Zc=100, Zf=50: t=0 %.4e, t=180 %.4e\n', bv2_0, bv2_f);
fprintf('Case 3 Zc=200, Zf=100: t=0 %.4e, t=180 %.4e\n', bv3_0, bv3_f);
fprintf('Case 2 vs case 1 at t=180: %.2f%%\n', 100 * (bv2_f - bv1_f) / max(bv1_f, eps));
fprintf('Case 3 vs case 1 at t=180: %.2f%%\n', 100 * (bv3_f - bv1_f) / max(bv1_f, eps));

fig_dir = fullfile(repo_root, 'output', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Figure 1: total biovolume vs depth at t=180
tot1 = sum(Y1_f, 2);
tot2 = sum(Y2_f, 2);
tot3 = sum(Y3_f, 2);
xmax = max([tot1(:); tot2(:); tot3(:)]);

f1 = figure;
subplot(1,3,1);
plot(tot1, z_plot, 'b-', 'LineWidth', 1.1);
set(gca, 'YDir', 'reverse');
xlim([0, xmax]);
xlabel('biovolume');
ylabel('depth (m)');
title('no grazing');

subplot(1,3,2);
plot(tot2, z_plot, 'r-', 'LineWidth', 1.1);
set(gca, 'YDir', 'reverse');
xlim([0, xmax]);
xlabel('biovolume');
ylabel('depth (m)');
title('Zc=100');

subplot(1,3,3);
plot(tot3, z_plot, 'k-', 'LineWidth', 1.1);
set(gca, 'YDir', 'reverse');
xlim([0, xmax]);
xlabel('biovolume');
ylabel('depth (m)');
title('Zc=200');

saveas(f1, fullfile(fig_dir, 'run_1d_zoo_production_depth.png'));

% Figure 2: size spectrum at t=180 for surface/mid/deep
sec = 1:cfg_base.n_sections;
ks = [1, 10, 20];

f2 = figure;
subplot(1,3,1);
semilogy(sec, Y1_f(ks(1),:), 'b-', sec, Y1_f(ks(2),:), 'r-', sec, Y1_f(ks(3),:), 'k-');
xlabel('section');
ylabel('biovolume');
title('no grazing');
legend({'surface','mid','deep'}, 'Location', 'best');

subplot(1,3,2);
semilogy(sec, Y2_f(ks(1),:), 'b-', sec, Y2_f(ks(2),:), 'r-', sec, Y2_f(ks(3),:), 'k-');
xlabel('section');
ylabel('biovolume');
title('Zc=100');
legend({'surface','mid','deep'}, 'Location', 'best');

subplot(1,3,3);
semilogy(sec, Y3_f(ks(1),:), 'b-', sec, Y3_f(ks(2),:), 'r-', sec, Y3_f(ks(3),:), 'k-');
xlabel('section');
ylabel('biovolume');
title('Zc=200');
legend({'surface','mid','deep'}, 'Location', 'best');

saveas(f2, fullfile(fig_dir, 'run_1d_zoo_production_spectrum.png'));
fprintf('Saved figures in %s\n', fig_dir);

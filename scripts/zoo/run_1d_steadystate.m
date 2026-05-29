% run_1d_steadystate
% 1-D steady-state check with and without zooplankton grazing.
% Production mode matches the 0-D slab: dY/dt = mu * Y(1,1), mu = 0.1 day^-1.

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z_plot   = col_grid.z_centers;

cfg_base = SimulationConfig( ...
    'n_sections', 20, ...
    't_final', 500, ...
    'delta_t', 1, ...
    'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law', ...
    'enable_coag', true, ...
    'enable_sinking', true, ...
    'enable_disagg', false, ...
    'proc_substeps', 20, ...
    'enable_surface_pp', true, ...
    'surface_pp_bin', 1, ...
    'surface_pp_mu', 0.1, ...
    'enable_zoo', true, ...
    'zoo_Zc', 100, ...
    'zoo_c', 1e-4, ...
    'zoo_Zf', 50, ...
    'zoo_s', 1e-4, ...
    'zoo_p', 0.3, ...
    'zoo_ic', 1);

% Case 1: no grazing
cfg1 = cfg_base.copy();
cfg1.enable_zoo = false;
sim1 = ColumnSimulation(cfg1, col_grid, profile);
out1 = sim1.run();

% Case 2: with grazing
cfg2 = cfg_base.copy();
cfg2.enable_zoo = true;
sim2 = ColumnSimulation(cfg2, col_grid, profile);
out2 = sim2.run();

Y1 = squeeze(out1.concentrations);
Y2 = squeeze(out2.concentrations);
t  = out1.time;
Y1f = squeeze(Y1(end, :, :));            % n_z x n_sec
Y2f = squeeze(Y2(end, :, :));            % n_z x n_sec

tot1 = squeeze(sum(sum(Y1, 2), 3));
tot2 = squeeze(sum(sum(Y2, 2), 3));

days_check = [0, 100, 200, 500];
idx = zeros(size(days_check));
for i = 1:numel(days_check)
    [~, idx(i)] = min(abs(t - days_check(i)));
end

fprintf('\nNo grazing totals\n');
fprintf('%-8s  %-14s\n', 'day', 'total_bv');
for i = 1:numel(days_check)
    fprintf('%-8d  %14.4e\n', days_check(i), tot1(idx(i)));
end

fprintf('\nWith grazing totals\n');
fprintf('%-8s  %-14s\n', 'day', 'total_bv');
for i = 1:numel(days_check)
    fprintf('%-8d  %14.4e\n', days_check(i), tot2(idx(i)));
end

delta500 = 100 * (tot2(idx(end)) - tot1(idx(end))) / max(tot1(idx(end)), eps);
fprintf('\nChange at day 500 (with vs no grazing): %.2f%%\n', delta500);

% Budget at day 500 for grazing case.
% Rate units: filter feeders [day^-1], flux feeders [day^-1] via w [m/day].
% ZooplanktonGrazing.graze() takes w in cm/s but converts back to m/day
% internally before computing rate_FL = w_mday * s * Zf.
% So the budget must use m/day directly.
w_mday = out2.w_z;                       % n_z x n_sec [m/day]
rate_ff  = cfg2.zoo_c * cfg2.zoo_Zc;     % filter feeder [day^-1], uniform
rate_fl  = w_mday * cfg2.zoo_s * cfg2.zoo_Zf;  % flux feeder [day^-1], size-dependent
rate_all = rate_ff + rate_fl;

gross_grazing = sum(rate_all .* Y2f, 'all');
fecal_return  = cfg2.zoo_p * gross_grazing;
net_grazing   = gross_grazing - fecal_return;

% Production rate — depends on mode.
use_mu = isprop(cfg2,'surface_pp_mu') && cfg2.surface_pp_mu > 0;
if use_mu
    prod_rate = cfg2.surface_pp_mu * Y2f(1, cfg2.surface_pp_bin);
else
    prod_rate = cfg2.surface_pp_rate;
end

% Complete bottom flux: advective + diffusive (uses ColumnTransport.bottomFluxDay).
% bottomFluxDay returns [m/day x bv] per bin. Divide by dz to get [bv/day]
% contribution to the column-integrated budget (d/dt sum Y_k).
flux_vec = ColumnTransport.bottomFluxDay(Y2f, w_mday, out2.profile.Kz, col_grid.dz);
bottom_flux_adv  = sum(w_mday(end, :) .* max(Y2f(end, :), 0)) / col_grid.dz;
bottom_flux_diff = (sum(flux_vec) - sum(w_mday(end, :) .* max(Y2f(end, :), 0))) / col_grid.dz;
bottom_flux_tot  = sum(flux_vec) / col_grid.dz;

% dBV/dt from time series (t=200 to t=500 window).
idx200 = find(t >= 200, 1);
dbv_dt = (tot2(end) - tot2(idx200)) / (t(end) - t(idx200));

fprintf('\nBudget at day 500 (with grazing)\n');
fprintf('production rate     : %.4e per day\n', prod_rate);
fprintf('gross grazing       : %.4e per day\n', gross_grazing);
fprintf('fecal return        : %.4e per day\n', fecal_return);
fprintf('net grazing         : %.4e per day\n', net_grazing);
fprintf('bottom flux (adv)   : %.4e per day\n', bottom_flux_adv);
fprintf('bottom flux (diff)  : %.4e per day\n', bottom_flux_diff);
fprintf('bottom flux (total) : %.4e per day\n', bottom_flux_tot);
fprintf('dBV/dt (t200-500)   : %.4e per day\n', dbv_dt);
fprintf('budget residual     : %.4e per day\n', prod_rate - net_grazing - bottom_flux_tot - dbv_dt);

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

% Figure 1: total biovolume vs time
f1 = figure;
plot(t, tot1, 'b-', 'LineWidth', 1.2); hold on;
plot(t, tot2, 'r-', 'LineWidth', 1.2); hold off;
xlabel('time (day)');
ylabel('total biovolume');
title('total bv vs time');
legend({'no grazing', 'with grazing'}, 'Location', 'best');
saveas(f1, fullfile(fig_dir, 'run_1d_steadystate_timeseries.png'));

% Figure 2: depth profile at day 500
prof1 = sum(Y1f, 2);
prof2 = sum(Y2f, 2);
xmax = max([prof1(:); prof2(:)]);

f2 = figure;
plot(prof1, z_plot, 'b-', 'LineWidth', 1.2); hold on;
plot(prof2, z_plot, 'r-', 'LineWidth', 1.2); hold off;
set(gca, 'YDir', 'reverse');
xlim([0, xmax]);
xlabel('biovolume');
ylabel('depth (m)');
title('depth profile t=500');
legend({'no grazing', 'with grazing'}, 'Location', 'best');
saveas(f2, fullfile(fig_dir, 'run_1d_steadystate_depth.png'));

fprintf('\nSaved:\n');
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_steadystate_timeseries.png'));
fprintf('  %s\n', fullfile(fig_dir, 'run_1d_steadystate_depth.png'));

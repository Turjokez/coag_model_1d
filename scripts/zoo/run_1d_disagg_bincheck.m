% run_1d_disagg_bincheck
% Check which bins have actual biovolume at t=180 (surface layer).
% Also check what D_max is relative to populated bins.

clear; close all; clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);

cfg = SimulationConfig( ...
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
    'surface_pp_mu', 0.1, ...
    'enable_zoo', true, ...
    'zoo_Zc', 100, ...
    'zoo_c', 1e-4, ...
    'zoo_Zf', 50, ...
    'zoo_s', 1e-4, ...
    'zoo_p', 0.3, ...
    'zoo_ic', 1);

sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();
Yf  = squeeze(out.concentrations(end, :, :));  % n_z x n_sec

% size grid: volume-equivalent diameter
g    = DerivedGrid(cfg);
r_v  = (0.75/pi * g.av_vol).^(1/3);
d_um = 2 * r_v * 1e4;  % cm -> um

% D_max at each layer
Dmax_A  = 9.39e-6;
eps_cm  = profile.eps;
eps_m   = eps_cm / 1e4;
dmax_mm = Dmax_A * eps_m.^(-1/4) * 1000;  % m -> mm

fprintf('Surface layer (z=25m) size distribution at t=180\n');
fprintf('%-4s  %-10s  %-14s  %-8s\n', 'bin', 'd_v (um)', 'biovolume', '% of total');
surf_tot = sum(Yf(1,:));
for b = 1:20
    pct = 100 * Yf(1,b) / max(surf_tot, eps);
    marker = '';
    if d_um(b) >= dmax_mm(1)*1000, marker = ' <- above D_max'; end
    fprintf('%-4d  %-10.1f  %-14.4e  %-8.2f%s\n', b, d_um(b), Yf(1,b), pct, marker);
end
fprintf('D_max at surface = %.2f mm\n', dmax_mm(1));
fprintf('Largest populated bin (bv > 1e-10): ');
populated = find(Yf(1,:) > 1e-10);
if ~isempty(populated)
    fprintf('bin %d  (d = %.1f um)\n', populated(end), d_um(populated(end)));
else
    fprintf('none\n');
end

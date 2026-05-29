% run_science_check.m
% Science check for the full 1-D integrated model.
%
% Four figures:
%   1. Depth profile over time (how particles fill the column)
%   2. Size spectrum at surface, mid, deep at t=365
%   3. D_max(z) vs grid ceiling (where disagg is active)
%   4. Budget: surface production, total bv, bottom flux over time

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

% --- setup ---
col_grid = ColumnGrid(1000, 20);
profile  = DepthProfile.typical(col_grid.z_centers);
z        = col_grid.z_centers;
n_z      = col_grid.n_z;
dz       = col_grid.dz;

cfg = SimulationConfig( ...
    'n_sections',        20, ...
    't_final',           365, ...
    'delta_t',           1, ...
    'sinking_law',       'kriest_8', ...
    'ds_kernel_mode',    'sinking_law', ...
    'enable_coag',       true, ...
    'enable_sinking',    true, ...
    'enable_disagg',     true, ...
    'disagg_mode',       'operator_split', ...
    'disagg_dmax_cm',    1.0, ...
    'proc_substeps',     20, ...
    'enable_surface_pp', true, ...
    'surface_pp_bin',    1, ...
    'surface_pp_mu',     0.1, ...
    'enable_zoo',        true, ...
    'zoo_Zc',            0.307, ...
    'zoo_Zf',            0.063, ...
    'zoo_c',             0.025, ...
    'zoo_s',             1.3e-5, ...
    'zoo_p',             0.5, ...
    'zoo_ic',            7);

fprintf('Running science check (t=365, n=20)...\n');
sim = ColumnSimulation(cfg, col_grid, profile);
out = sim.run();

Yhist  = out.concentrations;   % n_t x n_z x n_sec
t_out  = out.time;             % n_t x 1
w_z    = out.w_z;              % n_z x n_sec, m/day
n_t    = length(t_out);
n_sec  = cfg.n_sections;

% bin grid representative diameters (um)
d0    = 20;
d_k   = d0 * 2.^((2*(1:n_sec) - 1)/6);

% depth indices for spectrum figure
iz_surf = 1;                        % z = 25 m
iz_mid  = round(n_z / 2);          % z ~ 500 m
iz_deep = n_z;                      % z = 975 m

% --- derived quantities ---

% total biovolume per depth layer vs time: n_t x n_z
bv_layer = squeeze(sum(Yhist, 3));

% total column biovolume vs time
bv_total = sum(bv_layer, 2);

% surface production rate vs time: mu * phi_1 at layer 1, bin 1
mu      = cfg.surface_pp_mu;
pp_rate = mu * squeeze(Yhist(:, 1, 1));   % day^-1 * biovolume

% bottom flux vs time: sum_s w(n_z,s) * Y(t,n_z,s) [m/day * m^-3 = m^-2 day^-1]
% PP flux: mu * Y(1,1) * dz  [day^-1 * m^-3 * m = m^-2 day^-1]
% Both are now per unit horizontal area — directly comparable.
bflux = zeros(n_t, 1);
for ti = 1:n_t
    bflux(ti) = sum(w_z(end, :) .* squeeze(Yhist(ti, end, :))');
end
pp_flux = pp_rate * dz;   % area-integrated PP [m^-2 day^-1]

% D_max profile
Dmax_A   = 9.39e-6;
dmax_mm  = zeros(n_z, 1);
for k = 1:n_z
    eps_m     = profile.eps(k) / 1e4;
    dmax_mm(k) = Dmax_A * eps_m^(-0.25) * 1000;
end
ceil_n20_mm = d_k(n_sec) / 1000;   % top bin diameter in mm

% --- print summary ---
fprintf('\n--- Science check summary ---\n');
fprintf('  Surface layer (z=25m):\n');
fprintf('    t=365 bv     = %.4e\n', bv_layer(end, 1));
fprintf('    bv / total   = %.1f%%\n', 100*bv_layer(end,1)/bv_total(end));
fprintf('  Deep layer (z=975m):\n');
fprintf('    t=365 bv     = %.4e\n', bv_layer(end, end));
fprintf('    bv / total   = %.1f%%\n', 100*bv_layer(end,end)/bv_total(end));
fprintf('  Bottom flux at t=365 = %.4e bv m^-2 day^-1\n', bflux(end));
fprintf('  PP flux at t=365     = %.4e bv m^-2 day^-1\n', pp_flux(end));
fprintf('  flux / PP ratio      = %.3f  (~%.1f%% of production reaches 1000m)\n', ...
    bflux(end)/max(pp_flux(end),eps), 100*bflux(end)/max(pp_flux(end),eps));

% --- figures ---
fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---- Figure 1: depth profile at t = 1, 30, 90, 180, 365 ----
snap_days = [1, 30, 90, 180, 365];
colors    = {[0.7 0.7 0.7], [0 0.4 0.8], [0.8 0.2 0], [0.6 0 0.6], [0 0 0]};

f1 = figure;
for i = 1:length(snap_days)
    ti = min(snap_days(i) + 1, n_t);
    plot(bv_layer(ti, :), z, '-', 'Color', colors{i}, 'LineWidth', 1.2);
    hold on;
end
hold off;
set(gca, 'YDir', 'reverse');
xlabel('biovolume (m^{-3})');
ylabel('depth (m)');
legend({'t=1','t=30','t=90','t=180','t=365'}, 'Location', 'best');
title('depth profile over time');
saveas(f1, fullfile(fig_dir, 'sci_depth_profile.png'));

% ---- Figure 2: size spectrum at 3 depths at t=365 ----
Yfinal = squeeze(Yhist(end, :, :));   % n_z x n_sec

f2 = figure;
semilogy(d_k, Yfinal(iz_surf, :), 'b-', 'LineWidth', 1.2); hold on;
semilogy(d_k, Yfinal(iz_mid,  :), 'r-', 'LineWidth', 1.2);
semilogy(d_k, Yfinal(iz_deep, :), 'k-', 'LineWidth', 1.2);
hold off;
xlabel('diameter (\mum)');
ylabel('biovolume (m^{-3})');
legend({sprintf('z=%dm', z(iz_surf)), ...
        sprintf('z=%dm', z(iz_mid)),  ...
        sprintf('z=%dm', z(iz_deep))}, 'Location', 'best');
title('size spectrum at t=365');
saveas(f2, fullfile(fig_dir, 'sci_size_spectrum.png'));

% ---- Figure 3: D_max(z) with grid ceiling ----
f3 = figure;
plot(dmax_mm, z, 'k-', 'LineWidth', 1.5); hold on;
xline(ceil_n20_mm, 'r--', 'LineWidth', 1.2);
% shade region where disagg is inactive (D_max > ceiling)
iz_inactive = find(dmax_mm > ceil_n20_mm, 1, 'first');
if ~isempty(iz_inactive)
    y_patch = [z(iz_inactive); z(end); z(end); z(iz_inactive)];
    x_patch = [0; 0; max(dmax_mm)*1.1; max(dmax_mm)*1.1];
    patch(x_patch, y_patch, [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
plot(dmax_mm, z, 'k-', 'LineWidth', 1.5);   % redraw on top
hold off;
set(gca, 'YDir', 'reverse');
xlabel('size (mm)');
ylabel('depth (m)');
legend({'D_{max}(z)', 'n=20 ceiling (1.81 mm)', 'disagg inactive'}, 'Location', 'best');
title('disaggregation: D_{max} vs depth');
saveas(f3, fullfile(fig_dir, 'sci_dmax_physics.png'));

% ---- Figure 4: budget over time ----
% two panels: left = total bv + surface bv, right = production and bottom flux
f4 = figure;

subplot(1, 2, 1);
plot(t_out, bv_total,       'k-',  'LineWidth', 1.2); hold on;
plot(t_out, bv_layer(:, 1), 'b--', 'LineWidth', 1.0);
hold off;
xlabel('time (day)');
ylabel('biovolume (m^{-3})');
legend({'column total', 'surface layer'}, 'Location', 'best');
title('total and surface biovolume');

subplot(1, 2, 2);
plot(t_out, pp_flux, 'b-',  'LineWidth', 1.2); hold on;
plot(t_out, bflux,   'r--', 'LineWidth', 1.2);
hold off;
xlabel('time (day)');
ylabel('biovolume m^{-2} day^{-1}');
legend({'PP flux (surface)', 'bottom flux (1000m)'}, 'Location', 'best');
title('production vs sinking flux');

saveas(f4, fullfile(fig_dir, 'sci_budget.png'));

fprintf('\nFigures saved:\n');
fprintf('  sci_depth_profile.png\n');
fprintf('  sci_size_spectrum.png\n');
fprintf('  sci_dmax_physics.png\n');
fprintf('  sci_budget.png\n');

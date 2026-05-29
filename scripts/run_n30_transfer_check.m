% run_n30_transfer_check.m
% Run the full 1-D model at n=30, dt=0.4 day.
%
% Goal: check if deep disaggregation (now active with n=30 grid ceiling = 18 mm)
% brings transfer efficiency down from the n=20 result of 49% to a
% realistic ocean range (~1-15%).
%
% Transfer efficiency (TE) = bottom sinking flux / surface PP flux.
%
% n=20 known result: TE = 49% (deep disagg inactive below 225m).
% n=30 expected: D_max < grid ceiling throughout most of the column,
% so disagg breaks large aggregates in the deep --> less fast sinking.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

% --- grid and profile (same for both runs) ---
col_grid = ColumnGrid(1000, 20);   % 20 depth layers, dz = 50 m
profile  = DepthProfile.typical(col_grid.z_centers);
z        = col_grid.z_centers;
n_z      = col_grid.n_z;
dz       = col_grid.dz;

% --- base config (physics same as the reference n=20 run) ---
base_cfg = { ...
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
    'zoo_ic',            7};           % fecal to bin 8 (~115 um)

% n=20 known results (from run_science_check.m output, May 21 2026)
%   corrected TE = (bflux/dz) / (mu*Y11*dz) = 49%
te_n20 = 49.0;

% --- run n=30, dt=0.4 ---
fprintf('Running n=30, dt=0.4 day, t=365 days...\n');
fprintf('(n=30 grid ceiling = 18.2 mm, D_max active across full column)\n\n');

cfg30 = SimulationConfig(base_cfg{:}, ...
    'n_sections', 30, ...
    't_final',    365, ...
    'delta_t',    0.4);

sim30 = ColumnSimulation(cfg30, col_grid, profile);
out30 = sim30.run();

Yhist = out30.concentrations;   % n_t x n_z x n_sec
t_out = out30.time;
w_z   = out30.w_z;              % n_z x n_sec, m/day
n_t   = length(t_out);
mu    = cfg30.surface_pp_mu;

% --- transfer efficiency ---
% PP flux: density-dependent surface production, area-integrated [bv m^-2 day^-1]
pp_flux = mu * Yhist(end, 1, 1) * dz;

% bottom sinking flux: sum of w * Y at lowest layer, corrected to [bv m^-2 day^-1]
bflux_raw = sum(w_z(end,:) .* squeeze(Yhist(end, end, :))');
bflux     = bflux_raw / dz;

te_n30 = 100 * bflux / max(pp_flux, eps);

% --- total biovolume ---
bv_total = squeeze(sum(Yhist, [2 3]));

% --- max populated bin in deep layer ---
Yfinal  = squeeze(Yhist(end, :, :));   % n_z x n_sec
bv_deep = Yfinal(end, :);
maxbin_deep = find(bv_deep > 0, 1, 'last');

% --- D_max vs grid ceiling ---
Dmax_A    = 9.39e-6;
d_top_n30 = 20.0 * 2^((2*30 - 1)/6) / 1000;   % mm, grid ceiling for n=30
d_top_n20 = 20.0 * 2^((2*20 - 1)/6) / 1000;   % mm, grid ceiling for n=20

dmax_mm = zeros(n_z, 1);
for k = 1:n_z
    eps_m     = profile.eps(k) / 1e4;
    dmax_mm(k) = Dmax_A * eps_m^(-0.25) * 1000;
end

% depth where disagg becomes active for n=30
active30 = find(dmax_mm < d_top_n30);
if ~isempty(active30)
    disagg_active_to_n30 = z(active30(end));
else
    disagg_active_to_n30 = 0;
end

% --- print results ---
fprintf('=== Transfer Efficiency ===\n');
fprintf('  n=20  (dt=1.0 day):   TE = %.1f%%  (known)\n', te_n20);
fprintf('  n=30  (dt=0.4 day):   TE = %.1f%%\n', te_n30);
fprintf('  Ocean target range:   1 - 15%%\n\n');

fprintf('=== D_max vs grid ceiling ===\n');
fprintf('  n=20 ceiling = %.2f mm  (disagg active to ~225 m)\n', d_top_n20);
fprintf('  n=30 ceiling = %.2f mm  (disagg active to ~%.0f m)\n', d_top_n30, disagg_active_to_n30);
fprintf('  D_max at surface = %.2f mm\n', dmax_mm(1));
fprintf('  D_max at 975m    = %.2f mm\n\n', dmax_mm(end));

fprintf('=== n=30 run summary ===\n');
fprintf('  Total bv t=0:    %.4e\n', bv_total(1));
fprintf('  Total bv t=365:  %.4e\n', bv_total(end));
fprintf('  PP flux (t=365): %.4e bv m^-2 day^-1\n', pp_flux);
fprintf('  Bottom flux:     %.4e bv m^-2 day^-1\n', bflux);
fprintf('  Deep max bin:    %d / 30\n', maxbin_deep);

% --- figures ---
fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% bin diameters for x-axis
d_k30 = 20 * 2.^((2*(1:30) - 1)/6);   % um

% Figure 1: depth profile at t=365
f1 = figure;
bv_z = sum(Yfinal, 2);
plot(bv_z, z, 'k-', 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse');
xlabel('biovolume (m^{-3})');
ylabel('depth (m)');
title(sprintf('n=30 depth profile t=365  (TE = %.1f%%)', te_n30));
saveas(f1, fullfile(fig_dir, 'n30_depth_profile.png'));

% Figure 2: size spectrum at surface and deep, n=30
f2 = figure;
semilogy(d_k30, Yfinal(1, :),   'b-', 'LineWidth', 1.2); hold on;
semilogy(d_k30, Yfinal(end, :), 'k-', 'LineWidth', 1.2);
hold off;
xlabel('diameter (\mum)');
ylabel('biovolume');
legend({'surface (25m)', 'deep (975m)'}, 'Location', 'best');
title('n=30 size spectrum at t=365');
saveas(f2, fullfile(fig_dir, 'n30_size_spectrum.png'));

% Figure 3: D_max vs depth with n=20 and n=30 ceilings
f3 = figure;
plot(dmax_mm, z, 'k-', 'LineWidth', 1.5); hold on;
plot([d_top_n20 d_top_n20], [z(1) z(end)], 'b--', 'LineWidth', 1.2);
plot([d_top_n30 d_top_n30], [z(1) z(end)], 'r--', 'LineWidth', 1.2);
hold off;
set(gca, 'YDir', 'reverse');
xlabel('size (mm)');
ylabel('depth (m)');
legend({'D_{max}(z)', 'n=20 ceiling (1.81 mm)', 'n=30 ceiling (18.2 mm)'}, 'Location', 'best');
title('D_{max} vs grid ceiling');
saveas(f3, fullfile(fig_dir, 'n30_dmax_ceiling.png'));

fprintf('\nFigures saved to %s\n', fig_dir);
fprintf('  n30_depth_profile.png\n');
fprintf('  n30_size_spectrum.png\n');
fprintf('  n30_dmax_ceiling.png\n');

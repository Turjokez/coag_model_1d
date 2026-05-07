% run_may06_phase2_depth_scaling_test
% Test that depth-dependent kernel scaling changes aggregate size with depth.
%
% Two runs, same config:
%   flat  — DepthProfile with constant T, nu, eps at all depths
%   ocean — DepthProfile.typical() with realistic T/nu/eps gradients
%
% Checks:
%   1. No negative concentrations in either run.
%   2. Keep nonzero signal in column (avoid collapse-to-zero).
%   3. Depth effect: mean size differs between shallow and deep active layers.
%   4. Flat and ocean runs give different column PSD (scaling actually does something).
%
% Figure: PSD at surface / mid / deep for both runs at final time.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
if ~exist(log_dir, 'dir'), mkdir(log_dir); end

% ---- shared config ----
cfg = SimulationConfig();
cfg.n_sections     = 5;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.enable_coag    = true;
cfg.enable_disagg  = false;
cfg.enable_sinking = true;
cfg.enable_linear  = false;
cfg.r_to_rg        = 1.6;
cfg.alpha          = 0.003;  % weak stickiness: keeps coag-induced biovolume drift < ~2%
cfg.t_init         = 0;
cfg.t_final        = 30;
cfg.delta_t        = 1;
cfg.proc_substeps  = 20;     % substeps inside each delta_t for process-rate stability

% ---- depth grid ----
H   = 1000;
n_z = 20;
cgrid = ColumnGrid(H, n_z);

% ---- flat profile: constant T=20C, nu=0.01 cm2/s, eps=1e-4 cm2/s3 ----
z = cgrid.z_centers;
nz = numel(z);
prof_flat = DepthProfile( ...
    z, ...
    (20 + 273.15) .* ones(nz,1), ...   % T_K
    35             .* ones(nz,1), ...   % S
    1.025          .* ones(nz,1), ...   % rho [g/cm3]
    0.01           .* ones(nz,1), ...   % nu  [cm2/s]
    1e-4           .* ones(nz,1), ...   % eps [cm2/s3]
    1e-4           .* ones(nz,1));      % Kz  [m2/s]

% ---- ocean profile: realistic gradients ----
prof_ocean = DepthProfile.typical(z);
% typical() gives Kz in m2/s already; set a gentle diffusivity
prof_ocean.Kz = 1e-4 .* ones(nz,1);

% ---- run flat ----
fprintf('Running flat profile...\n');
sim_flat  = ColumnSimulation(cfg, cgrid, prof_flat);
res_flat  = sim_flat.run();

% ---- run ocean ----
fprintf('Running ocean profile...\n');
sim_ocean = ColumnSimulation(cfg, cgrid, prof_ocean);
res_ocean = sim_ocean.run();

Y_flat  = res_flat.concentrations;   % n_t x n_z x n_sec
Y_ocean = res_ocean.concentrations;
t       = res_flat.time;
n_t     = numel(t);
n_sec   = cfg.n_sections;
av_vol  = sim_flat.size_grid.av_vol(:)';  % 1 x n_sec

% ---- checks ----
fprintf('\n=== Phase 2 depth-scaling checks ===\n');

% 1. negatives
neg_flat  = sum(Y_flat(:)  < -1e-30);
neg_ocean = sum(Y_ocean(:) < -1e-30);
fprintf('neg_count flat:  %d   (expected 0)\n', neg_flat);
fprintf('neg_count ocean: %d   (expected 0)\n', neg_ocean);

% 2. keep nonzero signal in flat run (sanity)
bv0 = sum(squeeze(Y_flat(1,:,:))   .* repmat(av_vol, n_z, 1), 'all');
bvf = sum(squeeze(Y_flat(end,:,:)) .* repmat(av_vol, n_z, 1), 'all');
bv_change_pct = 100 * (bvf - bv0) / max(bv0, 1e-60);
fprintf('biovolume change flat: %.4f %%  (sectional drift; < ~2%% ok at alpha=0.003)\n', bv_change_pct);
fprintf('final/initial biovolume ratio flat: %.6f\n', bvf / max(bv0, 1e-60));

% 3. depth effect in ocean run: compare shallow and deepest active layer
col_sum_ocean = squeeze(sum(Y_ocean(end, :, :), 3));
active_idx = find(col_sum_ocean > 1e-20);
if numel(active_idx) < 2
    iz_shallow = 1;
    iz_deep_act = n_z;
else
    iz_shallow = active_idx(1);
    iz_deep_act = active_idx(end);
end
psd_shallow = squeeze(Y_ocean(end, iz_shallow, :))';  % 1 x n_sec
psd_deepact = squeeze(Y_ocean(end, iz_deep_act, :))';
mean_size_shallow = sum(psd_shallow .* av_vol) / max(sum(psd_shallow), 1e-60);
mean_size_deepact = sum(psd_deepact .* av_vol) / max(sum(psd_deepact), 1e-60);
fprintf('mean biovolume at shallow active layer (z=%.1f m): %.3e\n', z(iz_shallow), mean_size_shallow);
fprintf('mean biovolume at deep active layer (z=%.1f m):    %.3e\n', z(iz_deep_act), mean_size_deepact);

% 4. flat vs ocean differ at midwater
iz_mid = round(n_z / 2);
psd_flat_mid  = squeeze(Y_flat(end,  iz_mid, :))';
psd_ocean_mid = squeeze(Y_ocean(end, iz_mid, :))';
total_flat_mid  = sum(psd_flat_mid  .* av_vol);
total_ocean_mid = sum(psd_ocean_mid .* av_vol);
diff_pct = 100 * abs(total_flat_mid - total_ocean_mid) / max(total_flat_mid, 1e-60);
fprintf('midwater biovolume diff flat vs ocean: %.2f %%\n', diff_pct);

% ---- figure: PSD at surface / mid / deep, flat vs ocean ----
d_cm   = sim_flat.size_grid.getVolumeDiameters();   % bin diameters [cm]
d_mm   = d_cm .* 10;

iz_plot = [1, round(n_z/2), n_z];
depth_labels = {'surface', 'mid', 'deep'};
colors_flat  = [0.6 0.6 0.6; 0.4 0.4 0.4; 0.1 0.1 0.1];
colors_ocean = [0.2 0.6 1.0; 0.0 0.3 0.8; 0.0 0.1 0.5];

figure;
hold on;
for ip = 1:3
    iz = iz_plot(ip);
    pf = squeeze(Y_flat(end,  iz, :));
    po = squeeze(Y_ocean(end, iz, :));
    plot(d_mm, pf, '--', 'Color', colors_flat(ip,:),  'LineWidth', 1.4, ...
        'DisplayName', sprintf('flat %s',  depth_labels{ip}));
    plot(d_mm, po, '-',  'Color', colors_ocean(ip,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('ocean %s', depth_labels{ip}));
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Diameter [mm]');
ylabel('Concentration');
title('Phase 2: depth-scaled kernels — flat vs ocean profile');
legend('Location', 'southwest');

fname = fullfile(fig_dir, 'may06_phase2_depth_scaling.png');
saveas(gcf, fname);
fprintf('\nFigure saved: %s\n', fname);

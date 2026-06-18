% run_size_split_test.m
%
% Size-class breakdown: is the deep deficit uniform across sizes,
% or concentrated in large particles (500-2000 um)?
%
% Uses best config. Splits 100-2000 um into:
%   Small: 100-500 um
%   Large: 500-2000 um
%
% Focus on 275-475 m. (Ignore 125-175 m -- BC pile-up artifact.)
%
% Interpretation:
%   Both classes low at depth  -> magnitude / source problem
%   Only large class low       -> coagulation not growing particles fast enough
%                                 OR disagg breaking them too much

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Setup
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

k_bc   = 2;
k_plot = [4 6 8 10];          % 175, 275, 375, 475 m
z_mod  = col_grid.z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% model bin size masks
d_model_um = bc.d_model_um;
mask_mod_sm = d_model_um >= 100 & d_model_um <  500;
mask_mod_lg = d_model_um >= 500 & d_model_um < 2000;

% UVP size masks
mask_uvp_sm = uvpd.d_um >= 100 & uvpd.d_um <  500;
mask_uvp_lg = uvpd.d_um >= 500 & uvpd.d_um < 2000;
mask_z      = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ia, ib] = intersect(bc.dates, uvpd.dates);

% UVP reference by size class
phi_uvp_sm = zeros(numel(k_plot), 1);
phi_uvp_lg = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z) - z_mod(ki)));

        phi_u_sm = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp_sm));
        phi_u_lg = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp_lg));
        if size(phi_u_sm,1) < size(phi_u_sm,2), phi_u_sm = phi_u_sm'; end
        if size(phi_u_lg,1) < size(phi_u_lg,2), phi_u_lg = phi_u_lg'; end

        phi_uvp_sm(ki) = phi_uvp_sm(ki) + sum(phi_u_sm(iz, :));
        phi_uvp_lg(ki) = phi_uvp_lg(ki) + sum(phi_u_lg(iz, :));
    end
end
phi_uvp_sm = phi_uvp_sm / numel(ia);
phi_uvp_lg = phi_uvp_lg / numel(ia);

% ---------------------------------------------------------------
% 2. Run best config
% ---------------------------------------------------------------
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    if rel_change < spinup_tol
        fprintf('Converged at cycle %d\n', icyc);
        break;
    end
end

% final run
Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);
phi_mod_sm = zeros(numel(k_plot), 1);
phi_mod_lg = zeros(numel(k_plot), 1);
n_cast = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    for i_step = 1:steps_per_day
        Y(k_bc, :) = phi_bc_daily(i_day, :);
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        Y(k_bc, :) = phi_bc_daily(i_day, :);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Ytot = Y + Yfp;
        for ki = 1:numel(k_plot)
            phi_mod_sm(ki) = phi_mod_sm(ki) + sum(Ytot(k_plot(ki), mask_mod_sm));
            phi_mod_lg(ki) = phi_mod_lg(ki) + sum(Ytot(k_plot(ki), mask_mod_lg));
        end
        n_cast = n_cast + 1;
    end
end
phi_mod_sm = phi_mod_sm / max(n_cast, 1);
phi_mod_lg = phi_mod_lg / max(n_cast, 1);

ratio_sm = phi_mod_sm ./ max(phi_uvp_sm, 1e-20);
ratio_lg = phi_mod_lg ./ max(phi_uvp_lg, 1e-20);

% ---------------------------------------------------------------
% 3. Print
% ---------------------------------------------------------------
fprintf('\n--- BV ratio by size class ---\n');
fprintf('%-8s  %-14s  %-14s\n', 'Depth', '100-500 um', '500-2000 um');
for ki = 1:numel(k_plot)
    fprintf('%5.0f m    %5.2f           %5.2f\n', z_mod(ki), ratio_sm(ki), ratio_lg(ki));
end

% ---------------------------------------------------------------
% 4. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 12], 'Color', 'white');
hold on;
plot(ratio_sm, z_mod, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', '100-500 um');
plot(ratio_lg, z_mod, 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', '500-2000 um');
plot([1 1], [z_mod(1) z_mod(end)], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse', 'YLim', [150 510]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Size class split (best config)', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'size_split_ratio.png'));
fprintf('\nSaved size_split_ratio.png\n');

% ---------------------------------------------------------------
function cfg = cfg_best()
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.enable_zoo     = true;
cfg.enable_microbe = false;
cfg.enable_mining  = true;
cfg.alpha          = 0.10;
cfg.microbe_r0     = 0.0;
cfg.surface_pp_mu  = 0.0;
cfg.r_to_rg        = 1.6;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.mining_s       = 1.3e-5;
cfg.fp_alpha_cross = 0.5;
end

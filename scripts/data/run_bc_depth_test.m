% run_bc_depth_test.m
%
% BC sensitivity test: does forcing the boundary condition at a deeper
% depth improve the model/UVP ratio at 375-475 m?
%
% Logic:
%   If ratio at 475 m improves when BC is pushed to 200 m or 300 m,
%   then the surface->deep transport is the bottleneck (missing physics
%   above the BC depth). A 1-D fix is possible.
%
%   If ratio at 475 m does NOT improve even with BC at 300 m, then
%   the 1-D model cannot maintain the observed deep standing stock
%   even with perfect particle input at mid-depth. Structural limitation.
%
% Three cases:
%   BC at 100 m  (k_bc = 2)  -- current config
%   BC at 200 m  (k_bc = 4)
%   BC at 300 m  (k_bc = 6)
%
% Note: layers ABOVE the BC are unconstrained in each case and should
% be ignored. Only layers below the BC are the meaningful test.
%
% Base config: alpha=0.10, Da x5, r0=0, zoo on, mining on, microbe off.

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

% layer centers: 25, 75, 125, 175, 225, 275, 325, 375, 425, 475 ...
% BC at 100 m  -> k_bc=2 (center 75 m)
% BC at 200 m  -> k_bc=4 (center 175 m)
% BC at 300 m  -> k_bc=6 (center 275 m)
cases = { ...
    struct('bc_depth', 100, 'k_bc', 2, 'label', 'BC at 100 m (current)', 'color', 'k'), ...
    struct('bc_depth', 200, 'k_bc', 4, 'label', 'BC at 200 m',           'color', 'r'), ...
    struct('bc_depth', 300, 'k_bc', 6, 'label', 'BC at 300 m',           'color', 'b'), ...
};

k_plot = 2:10;
z_mod  = col_grid.z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_base();

% ---------------------------------------------------------------
% 2. Build UVP reference (same for all cases: mean over cast days)
%    Use the 100m case's uvpd struct; dates and depth grid are fixed.
% ---------------------------------------------------------------
bc_ref       = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
uvpd         = bc_ref.uvpd;
mask_uvp     = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z       = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ia, ib]  = intersect(bc_ref.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z) - z_mod(ki)));
        phi_u = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
        if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz, :));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);

% ---------------------------------------------------------------
% 3. Run each case
% ---------------------------------------------------------------
ratio_all = zeros(numel(k_plot), numel(cases));

for ic = 1:numel(cases)
    bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, ...
                       cases{ic}.bc_depth, k_plot);
    phi_bc_daily = bc.phi_bc_daily;
    n_days       = bc.n_days;
    k_bc         = cases{ic}.k_bc;

    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);

    % spinup
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
            fprintf('%s: converged at cycle %d\n', cases{ic}.label, icyc);
            break;
        end
    end

    % final run
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);
    phi_mod = zeros(numel(k_plot), 1);
    n_cast  = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            Ytot    = Y + Yfp;
            phi_mod = phi_mod + sum(Ytot(k_plot, :), 2);
            n_cast  = n_cast + 1;
        end
    end

    ratio_all(:, ic) = (phi_mod / max(n_cast, 1)) ./ phi_uvp_ref;
end

% ---------------------------------------------------------------
% 4. Print table
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('%-8s  %-18s  %-14s  %-14s\n', 'Depth', cases{1}.label, cases{2}.label, cases{3}.label);
for ki = 1:numel(k_plot)
    fprintf('%5.0f m    %5.2f               %5.2f           %5.2f\n', ...
        z_mod(ki), ratio_all(ki,1), ratio_all(ki,2), ratio_all(ki,3));
end

% ---------------------------------------------------------------
% 5. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 12], 'Color', 'white');
hold on;
for ic = 1:numel(cases)
    plot(ratio_all(:,ic), z_mod, [cases{ic}.color '-o'], ...
         'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end
% mark BC depths
bc_depths = [100 200 300];
for ic = 1:3
    plot([0 2.5], [bc_depths(ic) bc_depths(ic)], ...
         [cases{ic}.color '--'], 'LineWidth', 0.6, 'HandleVisibility', 'off');
end
plot([1 1], [z_mod(1) z_mod(end)], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse', 'YLim', [60 510], 'XLim', [0 2.5]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('BC depth sensitivity test', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'bc_depth_test.png'));
fprintf('\nSaved bc_depth_test.png\n');

% ---------------------------------------------------------------
function cfg = cfg_base()
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
cfg.enable_dvm     = false;
end

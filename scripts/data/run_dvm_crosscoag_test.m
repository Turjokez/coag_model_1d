% run_dvm_crosscoag_test.m
%
% Isolation test: does cross-coagulation remove deep DVM fecal before it
% can accumulate as UVP-visible standing stock?
%
% Cross-coag converts fecal pellets into marine snow (Y_fp -> Y) via
% differential settling. If fp_alpha_cross is too high, DVM fecal injected
% at depth gets absorbed into Y immediately and disappears from Y_fp.
%
% Three cases:
%   1. DVM off, cross-coag on  (baseline, best config)
%   2. DVM on,  cross-coag on  (alpha_cross = 0.5, previous test)
%   3. DVM on,  cross-coag off (alpha_cross = 0,   isolation test)
%
% If case 3 improves ratio at 375-475 m but case 2 does not, then
% cross-coag is absorbing the DVM source before it shows up in Y_fp.
%
% DVM params: p_dvm=1.0, f_fec=0.2 (strong rerouting, 80% goes deep).
% Base config: alpha=0.10, Da x5, r0=0, 100m BC, operator_split disagg.

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
k_plot = 2:10;
z_mod  = col_grid.z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

bc           = get_daily_bc_at_depth(uvp_file, cfg_make(false, 0.5), col_grid, 100, 3:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% UVP reference
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ia, ib] = intersect(bc.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    for ki = 1:numel(k_plot)
        z_target = z_mod(ki);
        [~, iz] = min(abs(uvpd.depth_m(mask_z) - z_target));
        phi_u = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
        if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz, :));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);

% ---------------------------------------------------------------
% 2. Three cases
% ---------------------------------------------------------------
cases = { ...
    struct('dvm', false, 'alpha_cross', 0.5, 'label', 'DVM off  / cross-coag on',  'color', 'k'), ...
    struct('dvm', true,  'alpha_cross', 0.5, 'label', 'DVM on   / cross-coag on',  'color', 'r'), ...
    struct('dvm', true,  'alpha_cross', 0.0, 'label', 'DVM on   / cross-coag off', 'color', 'b'), ...
};

ratio_all = zeros(numel(k_plot), 3);

for ic = 1:3
    cfg = cfg_make(cases{ic}.dvm, cases{ic}.alpha_cross);
    cfg.validate();

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
            Ytot = Y + Yfp;
            phi_mod = phi_mod + sum(Ytot(k_plot, :), 2);
            n_cast  = n_cast + 1;
        end
    end

    ratio_all(:, ic) = (phi_mod / max(n_cast,1)) ./ phi_uvp_ref;
end

% ---------------------------------------------------------------
% 3. Print table
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('Depth    DVM-off/CC-on   DVM-on/CC-on   DVM-on/CC-off\n');
for ki = 1:numel(k_plot)
    fprintf('%5.0f m   %5.2f           %5.2f          %5.2f\n', ...
        z_mod(ki), ratio_all(ki,1), ratio_all(ki,2), ratio_all(ki,3));
end

% ---------------------------------------------------------------
% 4. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 10], 'Color', 'white');
hold on;
for ic = 1:3
    plot(ratio_all(:,ic), z_mod, [cases{ic}.color '-o'], ...
         'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end
plot([1 1], [z_mod(1) z_mod(end)], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse', 'YLim', [60 510]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('DVM + cross-coag isolation test', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'dvm_crosscoag_test.png'));
fprintf('\nSaved dvm_crosscoag_test.png\n');

% ---------------------------------------------------------------
function cfg = cfg_make(enable_dvm, alpha_cross)
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;   % dummy for validate; actual Dmax comes from Dmax_A
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
cfg.fp_alpha_cross = alpha_cross;
% DVM: strong rerouting (80% of fecal goes deep)
cfg.enable_dvm     = enable_dvm;
cfg.dvm_p          = 1.0;
cfg.dvm_ffec       = 0.2;
cfg.dvm_feed_zmax  = 150;
cfg.dvm_zmin       = 300;
cfg.dvm_zmax       = 500;
end

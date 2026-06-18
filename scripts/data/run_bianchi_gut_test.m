% run_bianchi_gut_test.m
%
% Test gut-memory DVM: does a persistent gut pool with delayed release
% improve the deep model/UVP ratio compared to simple Archibald rerouting?
%
% Key difference from run_dvm_test.m:
%   Old: F_dvm resets to 0 every timestep (no memory)
%   New: G_gut persists between timesteps, clears with tau_gut = 0.25 day
%
% Three cases:
%   1. DVM off            -- best config baseline
%   2. DVM on, tau = 0.25 -- gut clears in ~6 hours
%   3. DVM on, tau = 1.0  -- gut clears in ~1 day (slower, spreads deeper)
%
% Strong rerouting params used: dvm_p=1.0, dvm_ffec=0.2 (frac_deep=0.80)
% Same as previous DVM tests for direct comparison.
%
% Check: BV ratio at 375-475 m, mass conservation.

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

bc           = get_daily_bc_at_depth(uvp_file, cfg_make(false, 0.25), col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% UVP reference: mean BV over cast days, 100-2000 um
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ia, ib] = intersect(bc.dates, uvpd.dates);

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
% 2. Three cases
% ---------------------------------------------------------------
cases = { ...
    struct('dvm', false, 'tau', 0.25, 'label', 'DVM off',        'color', 'k'), ...
    struct('dvm', true,  'tau', 0.25, 'label', 'Gut tau = 0.25 day (~6 h)', 'color', 'r'), ...
    struct('dvm', true,  'tau', 1.0,  'label', 'Gut tau = 1.0 day',  'color', 'b'), ...
};

ratio_all = zeros(numel(k_plot), 3);
total_bv  = zeros(3, 1);

for ic = 1:3
    cfg = cfg_make(cases{ic}.dvm, cases{ic}.tau);
    cfg.validate();

    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        % reset gut pool at start of each spinup cycle
        sim.rhs.G_gut = 0;
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
    Y         = zeros(col_grid.n_z, cfg.n_sections);
    Yfp       = zeros(col_grid.n_z, cfg.n_sections);
    sim.rhs.G_gut = 0;
    phi_mod   = zeros(numel(k_plot), 1);
    n_cast    = 0;

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
    total_bv(ic)     = mean(sum(Y + Yfp, 2));
end

% ---------------------------------------------------------------
% 3. Print results
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('%-8s  %-10s  %-22s  %-16s\n', 'Depth', cases{1}.label, cases{2}.label, cases{3}.label);
for ki = 1:numel(k_plot)
    fprintf('%5.0f m    %5.2f       %5.2f                  %5.2f\n', ...
        z_mod(ki), ratio_all(ki,1), ratio_all(ki,2), ratio_all(ki,3));
end

fprintf('\n--- Mass conservation ---\n');
for ic = 1:3
    chg = 100 * (total_bv(ic) - total_bv(1)) / total_bv(1);
    fprintf('%s: total BV = %.4e  (%.2f%% vs DVM off)\n', ...
        cases{ic}.label, total_bv(ic), chg);
end

% ---------------------------------------------------------------
% 4. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 12], 'Color', 'white');
hold on;
for ic = 1:3
    plot(ratio_all(:,ic), z_mod, [cases{ic}.color '-o'], ...
         'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end
plot([1 1], [z_mod(1) z_mod(end)], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse', 'YLim', [60 510]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Gut-memory DVM test', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'bianchi_gut_test.png'));
fprintf('\nSaved bianchi_gut_test.png\n');

% ---------------------------------------------------------------
function cfg = cfg_make(enable_dvm, tau_gut)
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
% DVM: strong rerouting (same as previous tests)
cfg.enable_dvm     = enable_dvm;
cfg.dvm_p          = 1.0;
cfg.dvm_ffec       = 0.2;    % frac_deep = 0.80
cfg.dvm_feed_zmax  = 150;
cfg.dvm_zmin       = 300;
cfg.dvm_zmax       = 500;
cfg.dvm_tau_gut    = tau_gut;
end

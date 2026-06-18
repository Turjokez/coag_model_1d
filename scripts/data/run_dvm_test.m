% run_dvm_test.m
%
% Test DVM fecal rerouting: compare enable_dvm=false vs true.
%
% DVM routes a fraction of surface-zone fecal production to deep Y_fp
% (300-500 m). Based on Archibald 2019. Parameters: p_dvm=0.5, f_fec=0.7.
%
% Three checks:
%   1. BV ratio (model/UVP) at each depth -- does deep ratio improve?
%   2. Size spectrum at 375 m -- does 200-1200 um tail improve?
%   3. Mass conservation -- DVM only reroutes, total BV should not increase.
%
% Note: DVM fecal goes to the same fecal bin as normal fecal (~115 um, bin 8).
% The deep size deficit (200-1200 um) may not fully fix. This is a first-pass
% source test, not a final size fix.
%
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

bc           = get_daily_bc_at_depth(uvp_file, cfg_base(false), col_grid, 100, 3:10);
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
        z_target = z_mod(ki);
        [~, iz] = min(abs(uvpd.depth_m(mask_z) - z_target));
        phi_u = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
        if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz, :));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);

% UVP spectrum at 375 m for size comparison (check 3)
k_375 = find(abs(z_mod - 375) == min(abs(z_mod - 375)), 1);
phi_uvp_375 = zeros(sum(mask_uvp), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    [~, iz] = min(abs(uvpd.depth_m(mask_z) - 375));
    phi_u = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
    if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
    phi_uvp_375 = phi_uvp_375 + phi_u(iz, :)';
end
phi_uvp_375 = phi_uvp_375 / numel(ia);

% ---------------------------------------------------------------
% 2. Run: DVM off and DVM on
% ---------------------------------------------------------------
dvm_flags  = {false, true};
labels     = {'DVM off', 'DVM on'};
colors     = {'k', 'r'};
ratio_all  = zeros(numel(k_plot), 2);
spec_375   = zeros(30, 2);   % model spectrum at 375 m
total_bv   = zeros(2, 1);    % total column BV for mass check

for ir = 1:2
    cfg = cfg_base(dvm_flags{ir});
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
            fprintf('%s: spinup converged at cycle %d\n', labels{ir}, icyc);
            break;
        end
    end

    % final run -- collect stats on cast days
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);
    phi_mod  = zeros(numel(k_plot), 1);
    spec_sum = zeros(cfg.n_sections, 1);
    n_cast   = 0;

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
            % spectrum at 375 m layer
            [~, k375_mod] = min(abs(col_grid.z_centers - 375));
            spec_sum = spec_sum + Ytot(k375_mod, :)';
            n_cast   = n_cast + 1;
        end
    end

    ratio_all(:, ir) = (phi_mod / max(n_cast,1)) ./ phi_uvp_ref;
    spec_375(:, ir)  = spec_sum / max(n_cast, 1);
    total_bv(ir)     = mean(sum(Y + Yfp, 2));
end

% ---------------------------------------------------------------
% 3. Print results
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('Depth    DVM-off   DVM-on\n');
for ki = 1:numel(k_plot)
    fprintf('%5.0f m   %5.2f     %5.2f\n', z_mod(ki), ratio_all(ki,1), ratio_all(ki,2));
end

fprintf('\n--- Mass conservation check ---\n');
fprintf('Total column BV: DVM-off = %.4e,  DVM-on = %.4e\n', total_bv(1), total_bv(2));
fprintf('Relative change: %.2f%%\n', 100*(total_bv(2)-total_bv(1))/total_bv(1));

% ---------------------------------------------------------------
% 4. Figure 1: BV ratio vs depth
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 9 10], 'Color', 'white');
hold on;
for ir = 1:2
    plot(ratio_all(:,ir), z_mod, [colors{ir} '-o'], ...
         'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', labels{ir});
end
plot([1 1], [z_mod(1) z_mod(end)], 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse', 'YLim', [60 510]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 7);
title('DVM test: BV ratio vs depth', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'dvm_bv_ratio.png'));
fprintf('\nSaved dvm_bv_ratio.png\n');

% ---------------------------------------------------------------
% 5. Figure 2: size spectrum at 375 m
% ---------------------------------------------------------------
cfg_tmp  = cfg_base(false);
sg       = DerivedGrid(cfg_tmp);
d_um_mod = sg.dcomb * 1e4;   % cm to um
mask_mod = d_um_mod >= 100 & d_um_mod < 2000;

d_uvp_375 = uvpd.d_um(mask_uvp);

figure('Units', 'centimeters', 'Position', [2 2 11 8], 'Color', 'white');
hold on;
plot(d_uvp_375, phi_uvp_375, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', 'UVP');
for ir = 1:2
    plot(d_um_mod(mask_mod), spec_375(mask_mod, ir), [colors{ir} '-s'], ...
         'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', labels{ir});
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('ESD (\mum)');
ylabel('\phi (ppmV)');
legend('Location', 'southwest', 'FontSize', 7);
title('375 m size spectrum', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'dvm_spectrum_375m.png'));
fprintf('Saved dvm_spectrum_375m.png\n');

% ---------------------------------------------------------------
function cfg = cfg_base(enable_dvm)
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
cfg.fp_alpha_cross = 0.5;
% DVM
cfg.enable_dvm     = enable_dvm;
cfg.dvm_p          = 1.0;
cfg.dvm_ffec       = 0.2;
cfg.dvm_feed_zmax  = 150;
cfg.dvm_zmin       = 300;
cfg.dvm_zmax       = 500;
end

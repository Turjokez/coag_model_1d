% run_fix_test.m
%
% Test fixes for two identified problems:
%   Problem 1: BC input too high (w * phi_UVP overestimates sinking flux)
%   Problem 2: Flux attenuates too fast with depth (Martin b = 1.72 vs 0.86)
%
% Configs tested:
%   1. Baseline  (bc_scale=1.0, zoo on,  disagg A*5)
%   2. BC fix    (bc_scale=0.2, zoo on,  disagg A*5)
%   3. Zoo off   (bc_scale=1.0, zoo off, disagg A*5)
%   4. Disagg A*1(bc_scale=1.0, zoo on,  disagg A*1)
%   5. BC + zoo  (bc_scale=0.2, zoo off, disagg A*5)
%   6. BC + Dagg (bc_scale=0.2, zoo on,  disagg A*1)
%   7. BC+zoo+Da (bc_scale=0.2, zoo off, disagg A*1)
%
% Metrics:
%   - Martin b (100-500 m, BV flux)     target: 0.858
%   - Durkin ratio at 125, 330, 500 m   target: 1.0
%
% Uses flux BC throughout.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path  = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file  = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
trap_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'sediment_trap_durkin', 'raw', ...
    'cb6a494508_EXPORTS_EXPORTSNA_JC214_classified_geltrap_particlefluxes.sb');
fig_dir   = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% Setup (shared across all configs)
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);
z_centers = col_grid.z_centers;

k_bc   = 2;
dz     = col_grid.dz;
n_z    = col_grid.n_z;
dt     = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg0 = cfg_base();
k_plot = 2:10;
bc   = get_daily_bc_at_depth(uvp_file, cfg0, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

d_cm  = bc.d_model_um * 1e-4;
d_m   = bc.d_model_um * 1e-6;
w_bin = (66 * d_cm .^ 0.62)';
V_bin = ((pi/6) * d_m .^ 3)';
mask_uvp = bc.d_model_um >= 100 & bc.d_model_um < 2000;

% Durkin trap reference
trap = load_durkin_flux(trap_file);
mask_trap   = trap.d_um >= 100 & trap.d_um < 2000;
trap_depths = [125, 330, 500];
id_trap     = zeros(1, 3);
for i = 1:3
    [~, id_trap(i)] = min(abs(trap.depths - trap_depths(i)));
end
z_trap     = trap.depths(id_trap);
trap_total = sum(trap.flux_agg(id_trap, mask_trap), 2, 'omitnan');  % 3x1

[~, k100] = min(abs(z_centers - 100));
[~, k500] = min(abs(z_centers - 500));
k_trap    = zeros(1, 3);
for i = 1:3
    [~, k_trap(i)] = min(abs(z_centers - z_trap(i)));
end

% ---------------------------------------------------------------
% Config list
% ---------------------------------------------------------------
labels = { ...
    'Baseline (bc=1.0, zoo, Da*5)', ...
    'BC fix   (bc=0.2, zoo, Da*5)', ...
    'Zoo off  (bc=1.0,   - , Da*5)', ...
    'Disagg*1 (bc=1.0, zoo, Da*1)', ...
    'BC+zoo   (bc=0.2,   - , Da*5)', ...
    'BC+Da    (bc=0.2, zoo, Da*1)', ...
    'BC+zoo+Da(bc=0.2,   - , Da*1)' };

bc_scales  = [1.0,  0.2,  1.0,  1.0,  0.2,  0.2,  0.2];
zoo_on     = [true, true, false,true, false,true, false];
disagg_A   = [5,    5,    5,    1,    5,    1,    1   ];

n_cfg = numel(labels);
b_model    = NaN(n_cfg, 1);
ratio_trap = NaN(n_cfg, 3);

% ---------------------------------------------------------------
% Run each config
% ---------------------------------------------------------------
for ic = 1:n_cfg
    fprintf('\n--- Config %d: %s ---\n', ic, labels{ic});

    cfg = cfg_base();
    cfg.enable_zoo     = zoo_on(ic);
    cfg.disagg_dmax_A  = 9.39e-6 * disagg_A(ic);
    phi_bc = phi_bc_daily * bc_scales(ic);

    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            flux_src = dt * (w_bin .* phi_bc(i_day, :)) / dz;
            for i_step = 1:steps_per_day
                Y(k_bc, :) = Y(k_bc, :) + flux_src;
                [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            end
        end
        phi_after = mean(sum(Y + Yfp, 2));
        if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
            fprintf('  Converged at cycle %d\n', icyc); break;
        end
    end

    % final run: accumulate BV flux and number flux on cast days
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    F_bv_sum = zeros(n_z, 1);
    F_n_sum  = zeros(n_z, 1);
    n_cast = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc, :) = Y(k_bc, :) + flux_src;
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            for k = 1:n_z
                F_bv_sum(k) = F_bv_sum(k) + sum(w_bin .* Y(k,:));
                F_n_sum(k)  = F_n_sum(k)  + sum(w_bin(mask_uvp) .* Y(k,mask_uvp) ./ V_bin(mask_uvp));
            end
            n_cast = n_cast + 1;
        end
    end
    F_bv = F_bv_sum / max(n_cast, 1);
    F_n  = F_n_sum  / max(n_cast, 1);

    % Martin b from BV flux
    Ft = F_bv(k100);  Fb = F_bv(k500);
    if Ft > 0 && Fb > 0
        b_model(ic) = -log(Fb/Ft) / log(z_centers(k500)/z_centers(k100));
    end

    % Durkin ratio at 3 depths (number flux, 100-2000 um)
    for i = 1:3
        ratio_trap(ic, i) = F_n(k_trap(i)) / max(trap_total(i), 1e-30);
    end

    fprintf('  Martin b = %.3f\n', b_model(ic));
    fprintf('  Durkin ratio: 125m=%.2f  330m=%.2f  500m=%.2f\n', ...
        ratio_trap(ic,1), ratio_trap(ic,2), ratio_trap(ic,3));
end

% ---------------------------------------------------------------
% Print summary table
% ---------------------------------------------------------------
fprintf('\n\n=== SUMMARY TABLE ===\n');
fprintf('%-32s  %-8s  %-8s  %-8s  %-8s\n', ...
    'Config', 'b_model', 'r_125m', 'r_330m', 'r_500m');
fprintf('%-32s  %-8s  %-8s  %-8s  %-8s\n', ...
    '', '(tgt:0.86)', '(tgt:1.0)', '(tgt:1.0)', '(tgt:1.0)');
for ic = 1:n_cfg
    fprintf('%-32s  %8.3f  %8.2f  %8.2f  %8.2f\n', ...
        labels{ic}, b_model(ic), ratio_trap(ic,1), ratio_trap(ic,2), ratio_trap(ic,3));
end

% ---------------------------------------------------------------
% Figure: Martin b and Durkin ratios
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 18 8], 'Color', 'white');

subplot(1, 4, 1);
barh(1:n_cfg, b_model, 'FaceColor', [0.7 0.7 0.7]);
hold on;
plot([0.858 0.858], [0 n_cfg+1], 'r--', 'LineWidth', 1.2);
set(gca, 'YTick', 1:n_cfg, 'YTickLabel', 1:n_cfg, 'FontSize', 6);
xlabel('Martin b');
title('b (tgt: 0.86)', 'FontWeight', 'normal');
xlim([0 2.5]);

for col = 1:3
    subplot(1, 4, col+1);
    barh(1:n_cfg, ratio_trap(:, col), 'FaceColor', [0.6 0.8 1.0]);
    hold on;
    plot([1 1], [0 n_cfg+1], 'r--', 'LineWidth', 1.2);
    set(gca, 'YTick', 1:n_cfg, 'YTickLabel', 1:n_cfg, 'FontSize', 6, 'XScale', 'log');
    xlabel('Model/Trap');
    title(sprintf('%dm (tgt:1)', z_trap(col)), 'FontWeight', 'normal');
end

% add config legend as text
annotation('textbox', [0.01 0.01 0.98 0.12], 'String', ...
    strjoin(arrayfun(@(i) sprintf('%d=%s',i,labels{i}), 1:n_cfg, 'UniformOutput',false), '  |  '), ...
    'FontSize', 4.5, 'EdgeColor', 'none', 'Interpreter', 'none');

saveas(gcf, fullfile(fig_dir, 'fix_sensitivity.png'));
fprintf('\nSaved fix_sensitivity.png\n');

% ---------------------------------------------------------------
function cfg = cfg_base()
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.disagg_dmax_A  = 9.39e-6 * 5;   % overridden per config
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.enable_zoo     = true;            % overridden per config
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

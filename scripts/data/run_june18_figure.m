% run_june18_figure.m
%
% Three-panel summary figure for June 18 report.
%
% Panel 1: Total biovolume vs depth — model vs UVP
% Panel 2: Model / UVP ratio vs depth
% Panel 3: BC fix effect — baseline vs bc*0.2 (flux attenuation profile)
%
% Shows the two problems:
%   (1) model/UVP ratio falls from ~0.8 to 0.13 with depth
%   (2) scaling down the BC fixes shallow but makes deep worse

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir  = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

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

cfg0   = cfg_base();
k_plot = 2:10;
bc     = get_daily_bc_at_depth(uvp_file, cfg0, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

d_cm  = bc.d_model_um * 1e-4;
w_bin = (66 * d_cm .^ 0.62)';
mask_uvp = bc.d_model_um >= 100 & bc.d_model_um < 2000;

z_mod = z_centers(k_plot);

% ---------------------------------------------------------------
% UVP reference biovolume at each depth
% ---------------------------------------------------------------
mask_uvp_raw = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z_uvp   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ia, ib]  = intersect(bc.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    phi_u  = squeeze(uvpd.phi(id_uvp, mask_z_uvp, mask_uvp_raw));
    if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z_uvp) - z_mod(ki)));
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz, :));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);
uvp_thresh  = 0.01 * max(phi_uvp_ref);
mask_ok     = phi_uvp_ref >= uvp_thresh;

% ---------------------------------------------------------------
% Run model: two configs — baseline and BC*0.2
% ---------------------------------------------------------------
bc_scales  = [1.0, 0.2];
cfg_labels = {'Baseline', 'BC \times0.2'};
line_styles = {'-', '--'};
phi_mod    = zeros(numel(k_plot), 2);
F_bv_full  = zeros(n_z, 2);

for ic = 1:2
    cfg  = cfg_base();
    phi_bc = phi_bc_daily * bc_scales(ic);
    sim  = ColumnSimulation(cfg, col_grid, prof);
    Y    = zeros(n_z, cfg.n_sections);
    Yfp  = zeros(n_z, cfg.n_sections);

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
            fprintf('Config %d converged at cycle %d\n', ic, icyc); break;
        end
    end

    % final run: accumulate on cast days
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    phi_sum  = zeros(numel(k_plot), 1);
    F_bv_sum = zeros(n_z, 1);
    n_cast   = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc, :) = Y(k_bc, :) + flux_src;
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            Ytot = Y + Yfp;
            phi_sum = phi_sum + sum(Ytot(k_plot, mask_uvp), 2);
            for k = 1:n_z
                F_bv_sum(k) = F_bv_sum(k) + sum(w_bin .* Y(k,:));
            end
            n_cast = n_cast + 1;
        end
    end
    phi_mod(:, ic)   = phi_sum   / max(n_cast, 1);
    F_bv_full(:, ic) = F_bv_sum  / max(n_cast, 1);
end

% ratio model/UVP (baseline)
ratio = NaN(numel(k_plot), 1);
ratio(mask_ok) = phi_mod(mask_ok, 1) ./ phi_uvp_ref(mask_ok);

% ---------------------------------------------------------------
% Figure
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 18 12], 'Color', 'white');

% --- Panel 1: total BV profile ---
subplot(1, 3, 1);
hold on;
plot(phi_uvp_ref, z_mod, 'b-o', 'MarkerSize', 4, 'LineWidth', 1.3, ...
    'DisplayName', 'UVP');
plot(phi_mod(:, 1), z_mod, 'k-', 'LineWidth', 1.3, ...
    'DisplayName', 'Model');
set(gca, 'YDir', 'reverse', 'XScale', 'log', ...
    'YLim', [60 510], 'FontSize', 7);
xlabel('Biovolume (m^3 m^{-3})');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Total BV (100–2000 \mum)', 'FontWeight', 'normal');
hold off;

% --- Panel 2: model/UVP ratio ---
subplot(1, 3, 2);
hold on;
plot(ratio, z_mod, 'k-o', 'MarkerSize', 4, 'LineWidth', 1.3);
plot([1 1], [60 510], 'r:', 'LineWidth', 1.0);
text(1.05, 100, 'ratio = 1', 'FontSize', 6, 'Color', 'r');
set(gca, 'YDir', 'reverse', 'XScale', 'log', ...
    'YLim', [60 510], 'XLim', [0.05 5], 'FontSize', 7);
xlabel('Model / UVP');
ylabel('Depth (m)');
title('Model / UVP ratio', 'FontWeight', 'normal');
hold off;

% --- Panel 3: BV flux profiles — baseline vs BC*0.2 ---
subplot(1, 3, 3);
hold on;
colors = {'k', [0.5 0.5 0.5]};
for ic = 1:2
    semilogy(F_bv_full(:, ic), z_centers, ...
        'Color', colors{ic}, 'LineStyle', line_styles{ic}, ...
        'LineWidth', 1.3, 'DisplayName', cfg_labels{ic});
end
set(gca, 'YDir', 'reverse', 'YLim', [60 550], 'FontSize', 7);
xlabel('BV flux (m^3 m^{-2} d^{-1})');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Flux profile — BC fix effect', 'FontWeight', 'normal');
hold off;

saveas(gcf, fullfile(fig_dir, 'june18_summary.png'));
fprintf('Saved june18_summary.png\n');

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
end

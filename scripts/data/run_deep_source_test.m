% run_deep_source_test.m
%
% Test C: Simple deep source diagnostic.
%
% Can a small sub-100m background production term close the deep gap?
%
% Source: S(z) = S0 * exp(-(z - 100) / H)  for z > 100 m
% Injected into small/mid bins (100-500 um) uniformly.
% Not meant to be realistic biology -- diagnostic only.
%
% Three S0 values tested: 1e-9, 5e-9, 1e-8
% H = 150 m (e-folding depth of sub-euphotic production)
%
% Uses flux BC throughout.
% Reports: deep ratio at 475 m, Martin b, total BV change.

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
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);

k_bc      = 2;
dz        = col_grid.dz;
n_z       = col_grid.n_z;
z_centers = col_grid.z_centers;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
k_plot_bc = 2:10;
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot_bc);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;
d_model_um   = bc.d_model_um;

d_cm  = d_model_um * 1e-4;
w_bin = (66 * d_cm .^ 0.62)';   % 1 x n_sec

% deep source: inject into 100-500 um bins
src_bins = d_model_um >= 100 & d_model_um < 500;
n_src    = sum(src_bins);

% source depth profile: exp decay below 100m
H_src    = 150;    % e-folding depth [m]
z0_src   = 100;
src_z    = zeros(n_z, 1);
for k = 1:n_z
    if z_centers(k) > z0_src
        src_z(k) = exp(-(z_centers(k) - z0_src) / H_src);
    end
end

% reference indices
[~, k100] = min(abs(z_centers - 100));
[~, k475] = min(abs(z_centers - 475));
[~, k500] = min(abs(z_centers - 500));

% UVP reference: 100-2000 um, mean over cast days
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
k_plot   = 2:10;
z_mod    = z_centers(k_plot);
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
uvp_thresh = 0.01 * max(phi_uvp_ref);
mask_ok    = phi_uvp_ref >= uvp_thresh;

% ---------------------------------------------------------------
% 2. Cases: flux BC only + three S0 values
% ---------------------------------------------------------------
S0_vals = [0, 1e-9, 5e-9, 1e-8];
labels  = {'No source (flux BC)', 'S0=1e-9', 'S0=5e-9', 'S0=1e-8'};
colors  = {'b', 'g', [1 0.5 0], 'r'};
n_cases = numel(S0_vals);

ratio_all    = NaN(numel(k_plot), n_cases);
flux_at_100  = zeros(n_cases, 1);
flux_at_500  = zeros(n_cases, 1);

for ic = 1:n_cases
    S0 = S0_vals(ic);
    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);

    % per-step source: S0 * src_z(k) / n_src, scaled by dt
    src_step = dt * S0 * src_z / n_src;   % n_z x 1

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            flux_src_bc = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
            for i_step = 1:steps_per_day
                % flux BC at surface
                Y(k_bc, :) = Y(k_bc, :) + flux_src_bc;
                % deep background source
                for k = 1:n_z
                    if src_step(k) > 0
                        Y(k, src_bins) = Y(k, src_bins) + src_step(k);
                    end
                end
                [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            end
        end
        phi_after  = mean(sum(Y + Yfp, 2));
        rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
        if rel_change < spinup_tol
            fprintf('%s: converged at cycle %d\n', labels{ic}, icyc);
            break;
        end
    end

    % final run
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    phi_mod      = zeros(numel(k_plot), 1);
    flux_acc_100 = 0;
    flux_acc_500 = 0;
    n_cast = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src_bc = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc, :) = Y(k_bc, :) + flux_src_bc;
            for k = 1:n_z
                if src_step(k) > 0
                    Y(k, src_bins) = Y(k, src_bins) + src_step(k);
                end
            end
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            Ytot = Y + Yfp;
            phi_mod      = phi_mod      + sum(Ytot(k_plot, :), 2);
            flux_acc_100 = flux_acc_100 + sum(w_bin .* Ytot(k100, :));
            flux_acc_500 = flux_acc_500 + sum(w_bin .* Ytot(k500, :));
            n_cast = n_cast + 1;
        end
    end

    phi_mod_avg     = phi_mod / max(n_cast, 1);
    flux_at_100(ic) = flux_acc_100 / max(n_cast, 1);
    flux_at_500(ic) = flux_acc_500 / max(n_cast, 1);

    r = NaN(numel(k_plot), 1);
    r(mask_ok) = phi_mod_avg(mask_ok) ./ phi_uvp_ref(mask_ok);
    ratio_all(:, ic) = r;
end

% ---------------------------------------------------------------
% 3. Print
% ---------------------------------------------------------------
fprintf('\n--- 475 m ratio and Martin b (100-500 m) ---\n');
fprintf('%-20s  %-10s  %-8s\n', 'Case', 'Ratio 475m', 'Martin b');
k475_idx = find(k_plot == k475);
for ic = 1:n_cases
    r475 = ratio_all(k475_idx, ic);
    F1   = flux_at_100(ic);
    F5   = flux_at_500(ic);
    if F1 > 0 && F5 > 0
        b = -log(F5/F1) / log(z_centers(k500)/z_centers(k100));
    else
        b = NaN;
    end
    r_str = sprintf('%.2f', r475); if isnan(r475), r_str = ' NaN'; end
    fprintf('%-20s  %s          %.3f\n', labels{ic}, r_str, b);
end

% ---------------------------------------------------------------
% 4. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 14], 'Color', 'white');
hold on;

patch([0 3 3 0], [350 350 510 510], [0.85 0.85 0.85], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');
plot([1 1], [z_mod(1) z_mod(end)], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

for ic = 1:n_cases
    c = colors{ic};
    plot(ratio_all(:, ic), z_mod, '-o', 'Color', c, ...
         'MarkerSize', 3, 'LineWidth', 1.1, 'DisplayName', labels{ic});
end

set(gca, 'YDir', 'reverse', 'YLim', [60 510], 'XLim', [0 2.5]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Deep source sensitivity', 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'deep_source_test.png'));
fprintf('\nSaved deep_source_test.png\n');

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

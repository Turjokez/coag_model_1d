% run_flux_bc_test.m
%
% Diagnostic: Dirichlet BC vs Flux BC.
%
% Current model uses a hard concentration reset every timestep:
%   Y(k_bc,:) = phi_bc    <- forces concentration, creates pile-up below
%
% Flux BC instead injects a source proportional to the sinking flux:
%   dY/dt|source = w_bin .* phi_bc / dz
%   Y(k_bc,:) += dt * w_bin .* phi_bc / dz    <- no hard reset, free evolution
%
% Physical interpretation:
%   - Dirichlet: "concentration at 75 m IS the UVP value, always"
%   - Flux BC:   "particles crossing 75 m downward = w * phi_uvp; column evolves freely"
%
% Questions:
%   1. Does flux BC reduce the 125-175 m hump?
%   2. Does deep ratio improve?
%   3. Is the overall profile shape more physical?
%
% Note: w_bin uses Kriest_8 (w = 66 * d_cm^0.62 m/day), same as model.

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
dz     = col_grid.dz;          % layer thickness [m]
k_plot = 2:10;
z_mod  = col_grid.z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;   % n_days x n_sec
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% sinking speeds for flux BC: Kriest_8, w = 66 * d_cm^0.62 [m/day]
d_cm  = bc.d_model_um * 1e-4;    % um -> cm
w_bin = (66 * d_cm .^ 0.62)';    % 1 x n_sec [m/day]

% UVP reference: 100-2000 um, mean over cast days
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

% UVP threshold for ratio (avoid divide-by-near-zero)
uvp_thresh = 0.01 * max(phi_uvp_ref);
mask_ok    = phi_uvp_ref >= uvp_thresh;

% ---------------------------------------------------------------
% 2. Two cases
% ---------------------------------------------------------------
cases = { ...
    struct('bc_type', 'dirichlet', 'label', 'Dirichlet BC (current)', 'color', 'k'), ...
    struct('bc_type', 'flux',      'label', 'Flux BC',                'color', 'b'), ...
};

ratio_all = NaN(numel(k_plot), 2);
total_bv  = zeros(2, 1);

for ic = 1:2
    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);

    use_flux = strcmp(cases{ic}.bc_type, 'flux');

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;  % source per step
            for i_step = 1:steps_per_day
                if use_flux
                    % flux BC: inject source, let concentration evolve
                    Y(k_bc, :) = Y(k_bc, :) + flux_src;
                    [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                else
                    % Dirichlet BC: hard reset before and after
                    Y(k_bc, :) = phi_bc_daily(i_day, :);
                    [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                    Y(k_bc, :) = phi_bc_daily(i_day, :);
                end
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
        flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            if use_flux
                Y(k_bc, :) = Y(k_bc, :) + flux_src;
                [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            else
                Y(k_bc, :) = phi_bc_daily(i_day, :);
                [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                Y(k_bc, :) = phi_bc_daily(i_day, :);
            end
        end
        if any(bc.dates(i_day) == uvpd.dates)
            Ytot    = Y + Yfp;
            phi_mod = phi_mod + sum(Ytot(k_plot, :), 2);
            n_cast  = n_cast + 1;
        end
    end

    phi_mod_avg = phi_mod / max(n_cast, 1);
    r = NaN(numel(k_plot), 1);
    r(mask_ok) = phi_mod_avg(mask_ok) ./ phi_uvp_ref(mask_ok);
    ratio_all(:, ic) = r;
    total_bv(ic) = mean(sum(Y + Yfp, 2));
end

% ---------------------------------------------------------------
% 3. Print
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('%-8s  %-28s  %-12s\n', 'Depth', cases{1}.label, cases{2}.label);
for ki = 1:numel(k_plot)
    r1 = ratio_all(ki, 1);
    r2 = ratio_all(ki, 2);
    s1 = sprintf('%5.2f', r1); if isnan(r1), s1 = '  NaN'; end
    s2 = sprintf('%5.2f', r2); if isnan(r2), s2 = '  NaN'; end
    fprintf('%5.0f m    %s                            %s\n', z_mod(ki), s1, s2);
end

fprintf('\n--- Total column BV ---\n');
for ic = 1:2
    chg = 100 * (total_bv(ic) - total_bv(1)) / total_bv(1);
    fprintf('%s: %.4e  (%+.1f%%)\n', cases{ic}.label, total_bv(ic), chg);
end

% ---------------------------------------------------------------
% 4. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 14], 'Color', 'white');
hold on;

% shade out-of-scope region
patch([0 3 3 0], [350 350 510 510], [0.85 0.85 0.85], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.4, 'HandleVisibility', 'off');

% ratio=1 line
plot([1 1], [z_mod(1) z_mod(end)], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

for ic = 1:2
    plot(ratio_all(:, ic), z_mod, [cases{ic}.color '-o'], ...
         'MarkerSize', 3.5, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end

set(gca, 'YDir', 'reverse', 'YLim', [60 510], 'XLim', [0 2.5]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Dirichlet vs Flux BC', 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'flux_bc_test.png'));
fprintf('\nSaved flux_bc_test.png\n');

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

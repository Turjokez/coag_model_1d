% run_bc_variability_test.m
%
% Does day-to-day variability in the surface BC matter for deep standing stock?
%
% The current model already uses a day-varying BC (actual UVP on cast days,
% nearest-neighbor fill on other days). This script compares:
%
%   Case 1: variable BC  -- current approach (day-by-day UVP)
%   Case 2: mean BC      -- constant spectrum = mean of all cast-day BCs
%
% If both give the same deep ratio -> temporal variability does not help.
% The deep deficit is structural. We cannot fix it by preserving episodic events.
%
% If variable BC gives higher deep ratio -> episodic flux events matter.
% Time-averaging the BC suppresses high-flux pulses that penetrate deeper.
%
% Focus: 275-475 m. (Ignore 125-175 m -- BC pile-up artifact.)

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

cfg = cfg_best();
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% variable BC: day-by-day (current)
phi_bc_var = bc.phi_bc_daily;

% mean BC: average only over actual cast days
[~, ia, ib] = intersect(bc.dates, uvpd.dates);
phi_cast_mean = mean(phi_bc_var(ia, :), 1);   % 1 x n_sec
phi_bc_mean   = repmat(phi_cast_mean, n_days, 1);

% UVP reference: 100-2000 um, mean over cast days
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z) - z_mod(ki)));
        phi_u = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
        if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz, :));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);

% also print BC variability info
phi_vis = sum(phi_bc_var(ia, :), 2);
fprintf('BC at 100m -- cast-day total BV:\n');
fprintf('  mean = %.4e,  std = %.4e,  max/mean = %.2f\n', ...
    mean(phi_vis), std(phi_vis), max(phi_vis)/mean(phi_vis));

% ---------------------------------------------------------------
% 2. Two cases
% ---------------------------------------------------------------
cases = { ...
    struct('phi', phi_bc_var,  'label', 'Variable BC (current)', 'color', 'k'), ...
    struct('phi', phi_bc_mean, 'label', 'Mean BC',               'color', 'r'), ...
};

ratio_all = zeros(numel(k_plot), 2);

for ic = 1:2
    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);
    phi_bc_use = cases{ic}.phi;

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            for i_step = 1:steps_per_day
                Y(k_bc, :) = phi_bc_use(i_day, :);
                [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                Y(k_bc, :) = phi_bc_use(i_day, :);
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
            Y(k_bc, :) = phi_bc_use(i_day, :);
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_use(i_day, :);
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
% 3. Print
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('%-8s  %-24s  %-12s\n', 'Depth', cases{1}.label, cases{2}.label);
for ki = 1:numel(k_plot)
    fprintf('%5.0f m    %5.2f                    %5.2f\n', ...
        z_mod(ki), ratio_all(ki,1), ratio_all(ki,2));
end

% ---------------------------------------------------------------
% 4. Plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 12], 'Color', 'white');
hold on;
for ic = 1:2
    plot(ratio_all(:,ic), z_mod, [cases{ic}.color '-o'], ...
         'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end
plot([1 1], [z_mod(1) z_mod(end)], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse', 'YLim', [60 510]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('BC variability test', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'bc_variability_test.png'));
fprintf('\nSaved bc_variability_test.png\n');

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

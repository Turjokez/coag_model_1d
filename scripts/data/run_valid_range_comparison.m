% run_valid_range_comparison.m
%
% Best-config model vs UVP comparison, with explicit valid-range marking.
%
% Valid range for 1-D comparison: 125-325 m
%   - 125-175 m: minor BC pile-up artifact (known, flagged)
%   - 175-325 m: cleanest comparison zone
%
% Out of scope: 375-475 m
%   - All six isolation tests show a structural gap (~3x under)
%   - Not explained by zoo, coag, DVM, BC variability, or size class
%   - Likely a non-local source outside 1-D model scope
%
% Figure: ratio profile with valid range shaded, out-of-scope region marked.
% Also prints mean ratio inside valid range (excluding 125-175 m artifact).

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
k_plot = 2:10;         % 125 to 475 m
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

% final run: accumulate on cast days
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
phi_mod = phi_mod / max(n_cast, 1);

% threshold: exclude depths where UVP reference is too sparse
% use 1% of maximum UVP reference as the floor
uvp_thresh = 0.01 * max(phi_uvp_ref);
mask_ok    = phi_uvp_ref >= uvp_thresh;

ratio      = NaN(numel(k_plot), 1);
ratio(mask_ok) = phi_mod(mask_ok) ./ phi_uvp_ref(mask_ok);

% ---------------------------------------------------------------
% 3. Print
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('%-8s  %-8s  %s\n', 'Depth', 'Ratio', 'Zone');
% k_plot = 2:10 -> z_mod = 75,125,175,225,275,325,375,425,475 m
zone_labels = {'BC artifact', 'BC artifact', 'valid', 'valid', 'valid', 'valid', ...
               'out of scope', 'out of scope', 'out of scope'};
for ki = 1:numel(k_plot)
    if isnan(ratio(ki))
        fprintf('%5.0f m    %5s     %s  [UVP sparse]\n', z_mod(ki), 'NaN', zone_labels{ki});
    else
        fprintf('%5.0f m    %5.2f     %s\n', z_mod(ki), ratio(ki), zone_labels{ki});
    end
end

% mean ratio in valid zone: 175-325 m, only where UVP is valid
% indices 3-6 in k_plot = 175, 225, 275, 325 m
valid_idx   = 3:6;
valid_mask  = mask_ok(valid_idx);
valid_ratio = ratio(valid_idx);
if any(valid_mask)
    fprintf('\nMean ratio 175-325 m (UVP-valid depths only): %.2f  (n=%d)\n', ...
        mean(valid_ratio(valid_mask)), sum(valid_mask));
else
    fprintf('\nMean ratio 175-325 m: no valid UVP depths in this range\n');
end

% ---------------------------------------------------------------
% 4. Plot
% ---------------------------------------------------------------
z_full = [z_mod(1); z_mod; z_mod(end)];   % for patch edges

figure('Units', 'centimeters', 'Position', [2 2 9 14], 'Color', 'white');
hold on;

% out-of-scope region (grey, bottom)
patch([0 3 3 0], [350 350 510 510], [0.85 0.85 0.85], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Out of 1-D scope');

% valid range (light green, middle): 175-325 m
patch([0 3 3 0], [150 150 350 350], [0.80 0.93 0.80], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'Valid range');

% BC artifact (light yellow, top)
patch([0 3 3 0], [100 100 200 200], [0.98 0.96 0.80], ...
      'EdgeColor', 'none', 'FaceAlpha', 0.5, 'DisplayName', 'BC artifact');

% ratio=1 line
plot([1 1], [100 510], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');

% ratio profile (NaN depths shown as open circles)
plot(ratio, z_mod, 'k-o', 'MarkerSize', 3.5, 'LineWidth', 1.3, 'DisplayName', 'Best config');
% mark sparse-UVP depths explicitly
plot(zeros(sum(~mask_ok),1), z_mod(~mask_ok), 'ko', 'MarkerSize', 4, ...
     'MarkerFaceColor', 'w', 'HandleVisibility', 'off');

set(gca, 'YDir', 'reverse', 'YLim', [100 510], 'XLim', [0 2.5]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('Valid range comparison', 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'valid_range_comparison.png'));
fprintf('\nSaved valid_range_comparison.png\n');

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

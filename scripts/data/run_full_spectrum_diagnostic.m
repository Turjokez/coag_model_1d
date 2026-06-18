% run_full_spectrum_diagnostic.m
%
% Full 30-bin spectrum at 375 m and 475 m.
% Tests H1: are particles piling up in >2000 um bins (invisible to UVP)?
%
% Shows:
%   - Model BV across ALL bins (not just 100-2000 um UVP range)
%   - UVP visible range (100-2000 um) overlaid as blue squares
%   - Grey shading marks UVP visible window
%   - Prints fraction of model BV in <100, 100-2000, >2000 um

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

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

d_model_um = bc.d_model_um;   % n_sec x 1, model bin centers [um]

% layer indices
k375 = 8;    % z_center = 375 m
k475 = 10;   % z_center = 475 m

% UVP size and depth masks
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_ok = uvpd.d_um(mask_uvp);
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
z_uvp    = uvpd.depth_m(mask_z);
[~, iz375] = min(abs(z_uvp - 375));
[~, iz475] = min(abs(z_uvp - 475));

[~, ia, ib] = intersect(bc.dates, uvpd.dates);

% ---------------------------------------------------------------
% 2. Pre-compute UVP reference at 375 m and 475 m
% ---------------------------------------------------------------
uvp375 = zeros(sum(mask_uvp), 1);
uvp475 = zeros(sum(mask_uvp), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    phi_u  = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
    if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
    uvp375 = uvp375 + phi_u(iz375, :)';
    uvp475 = uvp475 + phi_u(iz475, :)';
end
uvp375 = uvp375 / numel(ia);
uvp475 = uvp475 / numel(ia);

% ---------------------------------------------------------------
% 3. Spinup
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

% ---------------------------------------------------------------
% 4. Final run: accumulate model spectrum on cast days
% ---------------------------------------------------------------
Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);
spec375 = zeros(cfg.n_sections, 1);
spec475 = zeros(cfg.n_sections, 1);
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
        spec375 = spec375 + Ytot(k375, :)';
        spec475 = spec475 + Ytot(k475, :)';
        n_cast  = n_cast + 1;
    end
end
spec375 = spec375 / max(n_cast, 1);
spec475 = spec475 / max(n_cast, 1);

% ---------------------------------------------------------------
% 5. Print size fractions
% ---------------------------------------------------------------
mask_tiny = d_model_um <  100;
mask_mid  = d_model_um >= 100 & d_model_um < 2000;
mask_big  = d_model_um >= 2000;

fprintf('\n--- Model BV size fractions ---\n');
depth_labels = {'375 m', '475 m'};
specs_list   = {spec375, spec475};
for ip = 1:2
    sp  = specs_list{ip};
    tot = max(sum(sp), 1e-30);
    fprintf('%s:  <100 um = %.1f%%   100-2000 um = %.1f%%   >2000 um = %.1f%%\n', ...
        depth_labels{ip}, ...
        100*sum(sp(mask_tiny))/tot, ...
        100*sum(sp(mask_mid)) /tot, ...
        100*sum(sp(mask_big)) /tot);
end

% ---------------------------------------------------------------
% 6. Plot: 2-panel full spectrum
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 18 8], 'Color', 'white');

uvp_refs = {uvp375, uvp475};

for ip = 1:2
    subplot(1, 2, ip);
    hold on;

    % UVP visible window
    ylims = [1e-12 1e-4];
    patch([100 2000 2000 100], [ylims(1) ylims(1) ylims(2) ylims(2)], ...
          [0.9 0.9 0.9], 'EdgeColor', 'none', 'FaceAlpha', 0.4);

    % model: all 30 bins
    plot(d_model_um, specs_list{ip}, 'k-o', ...
         'MarkerSize', 2.5, 'LineWidth', 1.0, 'DisplayName', 'Model');

    % UVP: 100-2000 um only
    plot(d_uvp_ok, uvp_refs{ip}, 'b-s', ...
         'MarkerSize', 3, 'LineWidth', 1.0, 'DisplayName', 'UVP');

    set(gca, 'XScale', 'log', 'YScale', 'log', ...
        'XLim', [10 1e4], 'YLim', ylims, 'FontSize', 7);
    xlabel('Diameter (\mum)');
    ylabel('BV');
    title(depth_labels{ip}, 'FontWeight', 'normal');
    if ip == 1
        legend('Location', 'northeast', 'FontSize', 6);
    end
    hold off;
end

saveas(gcf, fullfile(fig_dir, 'full_spectrum_diagnostic.png'));
fprintf('\nSaved full_spectrum_diagnostic.png\n');

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

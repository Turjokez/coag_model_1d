% run_normalized_spectrum_test.m
%
% Test B: Normalized size spectrum shape comparison.
%
% Both model and UVP normalized to unit total BV at each depth.
% Removes absolute magnitude problem. Asks: even if amounts are wrong,
% is the SHAPE of the size distribution correct?
%
% If shape matches -> size evolution physics are right, only missing source.
% If shape mismatches -> model is producing wrong size distribution.
%
% Uses flux BC (more physical). Plots 3 depths: 125, 275, 475 m.

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
dz     = col_grid.dz;
n_z    = col_grid.n_z;

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
d_model_um   = bc.d_model_um;   % model bin centers [um]

d_cm  = d_model_um * 1e-4;
w_bin = (66 * d_cm .^ 0.62)';   % 1 x n_sec

% depths to compare: 125, 275, 475 m
z_compare = [125, 275, 475];
k_compare = zeros(1, 3);
for i = 1:3
    [~, k_compare(i)] = min(abs(col_grid.z_centers - z_compare(i)));
end

[~, ia, ib] = intersect(bc.dates, uvpd.dates);

% ---------------------------------------------------------------
% 2. UVP normalized spectra at comparison depths
% ---------------------------------------------------------------
mask_uvp_all = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_ok     = uvpd.d_um(mask_uvp_all);
mask_z_uvp   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
z_uvp        = uvpd.depth_m(mask_z_uvp);

uvp_spectra = zeros(numel(z_compare), sum(mask_uvp_all));
for m = 1:numel(ia)
    id_uvp = ib(m);
    phi_u  = squeeze(uvpd.phi(id_uvp, mask_z_uvp, mask_uvp_all));
    if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
    for i = 1:3
        [~, iz] = min(abs(z_uvp - z_compare(i)));
        uvp_spectra(i, :) = uvp_spectra(i, :) + phi_u(iz, :);
    end
end
uvp_spectra = uvp_spectra / numel(ia);

% normalize UVP spectra
uvp_norm = zeros(size(uvp_spectra));
for i = 1:3
    tot = sum(uvp_spectra(i, :));
    if tot > 0
        uvp_norm(i, :) = uvp_spectra(i, :) / tot;
    end
end

% ---------------------------------------------------------------
% 3. Run model (flux BC)
% ---------------------------------------------------------------
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc, :) = Y(k_bc, :) + flux_src;
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    if rel_change < spinup_tol
        fprintf('Converged at cycle %d\n', icyc);
        break;
    end
end

% final run: accumulate model spectra on cast days
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);
mod_spectra = zeros(3, cfg.n_sections);
n_cast = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc, :) = Y(k_bc, :) + flux_src;
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Ytot = Y + Yfp;
        for i = 1:3
            mod_spectra(i, :) = mod_spectra(i, :) + Ytot(k_compare(i), :);
        end
        n_cast = n_cast + 1;
    end
end
mod_spectra = mod_spectra / max(n_cast, 1);

% normalize model spectra (100-2000 um only for fair comparison)
mask_mod_uvp = d_model_um >= 100 & d_model_um < 2000;
mod_norm = zeros(3, cfg.n_sections);
for i = 1:3
    tot = sum(mod_spectra(i, mask_mod_uvp));
    if tot > 0
        mod_norm(i, mask_mod_uvp) = mod_spectra(i, mask_mod_uvp) / tot;
    end
end

% ---------------------------------------------------------------
% 4. Print: dominant bin at each depth
% ---------------------------------------------------------------
fprintf('\n--- Model dominant bin (100-2000 um) ---\n');
for i = 1:3
    sp = mod_norm(i, mask_mod_uvp);
    d  = d_model_um(mask_mod_uvp);
    [~, imax] = max(sp);
    fprintf('%5.0f m: peak at %.0f um  (%.1f%% of total)\n', ...
        z_compare(i), d(imax), 100*sp(imax)/sum(sp));
end

fprintf('\n--- UVP dominant bin (100-2000 um) ---\n');
for i = 1:3
    sp = uvp_norm(i, :);
    [~, imax] = max(sp);
    fprintf('%5.0f m: peak at %.0f um  (%.1f%% of total)\n', ...
        z_compare(i), d_uvp_ok(imax), 100*sp(imax)/sum(sp));
end

% ---------------------------------------------------------------
% 5. Plot: 3-panel normalized spectrum
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 18 7], 'Color', 'white');

for i = 1:3
    subplot(1, 3, i);
    hold on;

    % model: 100-2000 um bins
    plot(d_model_um(mask_mod_uvp), mod_norm(i, mask_mod_uvp), 'k-o', ...
         'MarkerSize', 2.5, 'LineWidth', 1.0, 'DisplayName', 'Model');

    % UVP
    plot(d_uvp_ok, uvp_norm(i, :), 'b-s', ...
         'MarkerSize', 3, 'LineWidth', 1.0, 'DisplayName', 'UVP');

    set(gca, 'XScale', 'log', 'YScale', 'log', ...
        'XLim', [90 2100], 'FontSize', 7);
    xlabel('Diameter (\mum)');
    if i == 1
        ylabel('Normalized BV');
        legend('Location', 'northeast', 'FontSize', 6);
    end
    title(sprintf('%d m', z_compare(i)), 'FontWeight', 'normal');
    hold off;
end

saveas(gcf, fullfile(fig_dir, 'normalized_spectrum_test.png'));
fprintf('\nSaved normalized_spectrum_test.png\n');

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

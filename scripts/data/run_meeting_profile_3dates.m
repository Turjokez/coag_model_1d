% run_meeting_profile_3dates.m
%
% Clean meeting figure: BV vs depth for 3 cast dates.
% UVP (black dots) vs model inverse-fit (red line).
% 3 side-by-side panels, same axes, no colorbar.
%
% Parameters: alpha=0.093, bc_scale=0.420, r0=0.014 (100m BC fit)
% Dates: 05-15, 05-21, 05-27  (early / mid / late cruise)
%
% Saves: docs/figures/meeting_profile_3dates.png

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Config
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

k_bc   = 2;
k_plot = 2:10;    % z centers 75-475 m
z_mod  = col_grid.z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1/dt);
spinup_tol    = 0.01;
max_cycles    = 80;
n_z           = col_grid.n_z;
dz            = col_grid.dz;

bc_scale = 0.420;

cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.enable_zoo     = true;
cfg.enable_microbe = true;
cfg.enable_mining  = true;
cfg.alpha          = 0.093;
cfg.microbe_r0     = 0.014;
cfg.surface_pp_mu  = 0.0;
cfg.r_to_rg        = 1.6;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.mining_s       = 1.3e-5;
cfg.fp_alpha_cross = 0.5;

% ---------------------------------------------------------------
% 2. BC + bin setup
% ---------------------------------------------------------------
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp  = ColumnSimulation(cfg, col_grid, prof);
n_sec    = cfg.n_sections;
d_cm     = sim_tmp.size_grid.dcomb(:)';
w_bin    = 66 * d_cm .^ 0.62;
d_um_mod = d_cm * 1e4;
bin_mask = d_um_mod >= 100 & d_um_mod < 2000;

% UVP bin mask 100-2000 um
uvp_bin_mask = uvpd.d_um >= 100 & uvpd.d_um < 2000;

% ---------------------------------------------------------------
% 3. Spinup
% ---------------------------------------------------------------
fprintf('Running spinup...\n');
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) * bc_scale / dz;
        for i_step = 1:steps_per_day
            Y(k_bc,:) = Y(k_bc,:) + flux_src;
            [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
        end
    end
    phi_after = mean(sum(Y + Yfp, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('Converged at cycle %d\n', icyc); break;
    end
end

% ---------------------------------------------------------------
% 4. Final run: continue from spun-up state, save daily BV at k_plot depths
% ---------------------------------------------------------------
% do NOT reset Y here — keep the spun-up equilibrium
bv_daily = zeros(numel(k_plot), n_days);   % total BV 100-2000 um

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) * bc_scale / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
    end
    Ytot = Y(k_plot,:) + Yfp(k_plot,:);
    bv_daily(:, i_day) = sum(Ytot(:, bin_mask), 2);
end
fprintf('Done\n');

% ---------------------------------------------------------------
% 5. Match dates + extract UVP BV per cast
% ---------------------------------------------------------------
sel_dates = [20210515, 20210521, 20210527];
date_lbls = {'05-15', '05-21', '05-27'};

[~, ia, ib] = intersect(bc.dates, uvpd.dates);
cast_dates = bc.dates(ia);

% UVP: total BV 100-2000 um per depth, per cast
n_uvp_z  = size(uvpd.phi, 2);
bv_uvp   = squeeze(sum(uvpd.phi(:, :, uvp_bin_mask), 3));   % [n_cast x n_uvp_z]
z_uvp    = uvpd.depth_m;
mask_z   = z_uvp >= 75 & z_uvp <= 510;
z_uvp_pl = z_uvp(mask_z);

% ---------------------------------------------------------------
% 6. Figure: 1 row x 3 panels
% ---------------------------------------------------------------
fs = 7;
xl = [3e-7 3e-5];   % BV range [m^3 m^-3]

figure('Units','centimeters','Position',[2 2 14 9],'Color','white');

for m = 1:3
    target_date = sel_dates(m);

    % model: find day index
    [~, ic] = intersect(cast_dates, target_date);
    if isempty(ic)
        fprintf('Date %d not in cast days, skipping\n', target_date); continue;
    end
    id_mod = ia(ic);
    bv_mod_day = bv_daily(:, id_mod);

    % UVP: find date index
    id_uvp = ib(ic);
    bv_uvp_day = bv_uvp(id_uvp, mask_z);

    subplot(1, 3, m);
    hold on;
    plot(bv_uvp_day, z_uvp_pl, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'k', ...
        'DisplayName', 'UVP');
    plot(bv_mod_day, z_mod, 'r-', 'LineWidth', 1.4, ...
        'DisplayName', 'model');
    set(gca, 'YDir', 'reverse', 'XScale', 'log', 'FontSize', fs, 'Box', 'on', ...
        'YLim', [75 510], 'XLim', xl);
    xlabel('BV 100-2000 \mum (m^3 m^{-3})', 'FontSize', fs);
    if m == 1
        ylabel('depth (m)', 'FontSize', fs);
        legend('Location', 'southeast', 'FontSize', fs, 'Box', 'off');
    else
        set(gca, 'YTickLabel', {});
    end
    title(date_lbls{m}, 'FontWeight', 'normal', 'FontSize', fs);
end

saveas(gcf, fullfile(fig_dir, 'meeting_profile_3dates.png'));
fprintf('Saved meeting_profile_3dates.png\n');

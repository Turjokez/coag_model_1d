% run_report_fig2_spectrum.m
%
% Figure 2 for combined report.
% Size spectrum (BV per bin) at 3 depths: 125 m, 325 m, 475 m.
% Model vs UVP on same axes.
%
% This shows WHERE in the size distribution the model is wrong:
%   - If the shape matches but the level is off -> total mass problem
%   - If the shape is different -> size bias problem
%
% Saves: report_fig2_spectrum.png

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir  = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- setup ---
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);
z_centers = col_grid.z_centers;
n_z       = col_grid.n_z;
dz        = col_grid.dz;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 80;
k_plot        = 2:10;
k_bc          = 2;

% depths to show: 125, 325, 475 m
plot_depths = [125, 325, 475];
n_pd        = numel(plot_depths);
k_pd        = zeros(1, n_pd);
for i = 1:n_pd
    [~, k_pd(i)] = min(abs(z_centers - plot_depths(i)));
end

% --- load UVP ---
cfg0 = cfg_base();
bc   = get_daily_bc_at_depth(uvp_file, cfg0, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp = ColumnSimulation(cfg0, col_grid, prof);
d_cm    = sim_tmp.size_grid.dcomb(:)';
d_um    = d_cm * 1e4;
w_bin   = (66 * d_cm .^ 0.62);

% model bin diameters (for x-axis)
d_mod_um = d_um;

% UVP bin diameters
d_uvp_um = uvpd.d_um;
mask_uvp_raw = d_uvp_um >= 100 & d_uvp_um < 2000;
d_uvp_plot   = d_uvp_um(mask_uvp_raw);

% UVP spectrum at each target depth, cast-day average
mask_z_uvp  = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ~, ib]  = intersect(bc.dates, uvpd.dates);

phi_uvp_spec = zeros(n_pd, sum(mask_uvp_raw));  % depth x bin
for m = 1:numel(ib)
    phi_u = squeeze(uvpd.phi(ib(m), mask_z_uvp, mask_uvp_raw));
    if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
    for i = 1:n_pd
        [~, iz] = min(abs(uvpd.depth_m(mask_z_uvp) - plot_depths(i)));
        phi_uvp_spec(i,:) = phi_uvp_spec(i,:) + phi_u(iz,:);
    end
end
phi_uvp_spec = phi_uvp_spec / max(numel(ib), 1);

% --- run model ---
fprintf('Running model...\n');
cfg = cfg_base();
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
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

% accumulate model spectrum on cast days
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);
phi_mod_spec = zeros(n_pd, cfg.n_sections);
n_cast = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Ytot = Y + Yfp;
        for i = 1:n_pd
            phi_mod_spec(i,:) = phi_mod_spec(i,:) + Ytot(k_pd(i),:);
        end
        n_cast = n_cast + 1;
    end
end
phi_mod_spec = phi_mod_spec / max(n_cast, 1);

% find shared y-limits across all panels
all_vals = [phi_uvp_spec(:); phi_mod_spec(:)];
all_vals = all_vals(all_vals > 0);
ylim_lo  = 10^floor(log10(min(all_vals)));
ylim_hi  = 10^ceil(log10(max(all_vals)));

% --- figure: 3 panels ---
fs = 7;
figure('Units','centimeters','Position',[2 2 18 7],'Color','white');

for i = 1:n_pd
    subplot(1, n_pd, i);
    hold on;
    plot(d_uvp_plot,   phi_uvp_spec(i,:), 'b-o', ...
        'MarkerSize',3,'LineWidth',1.2);
    plot(d_mod_um,     phi_mod_spec(i,:), 'k-',  ...
        'LineWidth',1.2);
    set(gca,'XScale','log','YScale','log', ...
        'XLim',[100 2000],'YLim',[ylim_lo ylim_hi], ...
        'FontSize',fs,'Box','off');
    xlabel('ESD (\mum)','FontSize',fs);
    if i == 1
        ylabel('BV (m^3 m^{-3})','FontSize',fs);
        legend({'UVP','Model'},'Location','northeast','FontSize',fs,'Box','off');
    end
    title(sprintf('%d m', plot_depths(i)),'FontWeight','normal','FontSize',fs);
    hold off;
end

saveas(gcf, fullfile(fig_dir, 'report_fig2_spectrum.png'));
fprintf('Saved report_fig2_spectrum.png\n');

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

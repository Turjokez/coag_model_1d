% run_cast_spectrum_inverse_fit.m
%
% Same cast-by-cast 2D spectrum figure as run_report_fig4_cast2d.m,
% but using the inverse-fit parameters (100m BC):
%   alpha = 0.093, bc_scale = 0.420, r0 = 0.014 day^-1
%
% Saves: docs/figures/cast_spectrum_inverse_fit.png

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
% 1. Config + grid
% ---------------------------------------------------------------
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);
z_centers = col_grid.z_centers;

k_bc   = 2;
k_plot = 2:10;     % model layers 2-10 (centers 75-475 m)
z_mod  = z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 80;
n_z           = col_grid.n_z;
dz            = col_grid.dz;

% inverse-fit parameters
bc_scale = 0.420;

cfg = SimulationConfig();
cfg.n_sections    = 30;
cfg.sinking_law   = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.disagg_mode   = 'operator_split';
cfg.disagg_dmax_A = 9.39e-6 * 5;
cfg.enable_coag   = true;
cfg.enable_disagg = true;
cfg.enable_zoo    = true;
cfg.enable_microbe = true;
cfg.enable_mining = true;
cfg.alpha         = 0.093;
cfg.microbe_r0    = 0.014;
cfg.surface_pp_mu = 0.0;
cfg.r_to_rg       = 1.6;
cfg.zoo_c         = 0.025;
cfg.zoo_s         = 1.3e-5;
cfg.zoo_p         = 0.5;
cfg.zoo_ic        = 7;
cfg.mining_s      = 1.3e-5;
cfg.fp_alpha_cross = 0.5;

% ---------------------------------------------------------------
% 2. BC + UVP
% ---------------------------------------------------------------
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;   % [n_days x n_sec]
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp = ColumnSimulation(cfg, col_grid, prof);
n_sec   = cfg.n_sections;
d_cm    = sim_tmp.size_grid.dcomb(:)';
w_bin   = 66 * d_cm .^ 0.62;

% model bin edges [mm]
d_mod_mm = d_cm * 10;
d_edges  = zeros(n_sec + 1, 1);
d_edges(1)       = d_mod_mm(1)^2 / d_mod_mm(2);
d_edges(n_sec+1) = d_mod_mm(n_sec)^2 / d_mod_mm(n_sec-1);
for k = 2:n_sec
    d_edges(k) = sqrt(d_mod_mm(k-1) * d_mod_mm(k));
end
dw_mod_mm = diff(d_edges);

% UVP bins 100-10000 um (0.1-10 mm, matching Siegel Fig 2a range)
mask_uvp  = uvpd.d_um >= 100 & uvpd.d_um <= 10000;
d_uvp_mm  = uvpd.d_um(mask_uvp) / 1000;
dw_uvp_mm = uvpd.dw(mask_uvp)   / 1000;

% UVP depth: 0-500 m
mask_z_uvp = uvpd.depth_m >= 0 & uvpd.depth_m <= 510;
z_uvp      = uvpd.depth_m(mask_z_uvp);

% remap model bins to UVP bins
n_uvp = sum(mask_uvp);
remap_frac = zeros(n_uvp, n_sec);
for j = 1:n_uvp
    uvp_lo = d_uvp_mm(j) - dw_uvp_mm(j);
    uvp_hi = d_uvp_mm(j);
    for k = 1:n_sec
        lo = max(d_edges(k),   uvp_lo);
        hi = min(d_edges(k+1), uvp_hi);
        if hi > lo
            remap_frac(j, k) = (hi - lo) / dw_mod_mm(k);
        end
    end
end

% ---------------------------------------------------------------
% 3. Spinup (flux BC with bc_scale)
% ---------------------------------------------------------------
fprintf('Running spinup (alpha=%.3f, bc_scale=%.3f, r0=%.4f)...\n', ...
    cfg.alpha, bc_scale, cfg.microbe_r0);
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
% 4. Final run: save daily snapshots
% ---------------------------------------------------------------
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);
Y_daily = zeros(numel(k_plot), n_sec, n_days);

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) * bc_scale / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
    end
    Y_daily(:,:,i_day) = Y(k_plot,:) + Yfp(k_plot,:);
end
fprintf('Model run complete\n');

% ---------------------------------------------------------------
% 5. Select 9 on-station dates (same as Siegel Fig 2a)
% ---------------------------------------------------------------
sel_dates = [20210514, 20210515, 20210516, 20210519, ...
             20210521, 20210522, 20210525, 20210527, 20210529];

[~, ia, ib] = intersect(bc.dates, uvpd.dates);
cast_dates   = bc.dates(ia);
uvp_cast_idx = ib;

[~, isel_c, ~] = intersect(cast_dates, sel_dates);
n_sel = numel(isel_c);
fprintf('Selected %d dates\n', n_sel);

date_lbls = cell(n_sel, 1);
for m = 1:n_sel
    ds = num2str(cast_dates(isel_c(m)));
    date_lbls{m} = [ds(5:6) '-' ds(7:8)];
end

% ---------------------------------------------------------------
% 6. Figure
% ---------------------------------------------------------------
cmin = -1;   cmax = 1;
xlim_log  = [log10(0.09) log10(12)];
xtick_pos = log10([0.1 1.0 10]);
xtick_lbl = {'0.1', '1', '10'};

x_log = log10(d_uvp_mm);

fig_w = max(n_sel * 2.0 + 1.5, 16);
figure('Units','centimeters','Position',[1 1 fig_w 10],'Color','white');
colormap(jet);

for m = 1:n_sel
    id_cast = isel_c(m);
    id_mod  = ia(id_cast);
    id_uvp  = uvp_cast_idx(id_cast);

    % UVP spectrum [ppmV/mm], 0-500 m
    phi_all = squeeze(uvpd.phi(id_uvp, :, :));
    phi_u   = phi_all(mask_z_uvp, mask_uvp);
    phi_u(isnan(phi_u)) = 0;
    S_u = bsxfun(@rdivide, phi_u, dw_uvp_mm) * 1e6;
    S_u(S_u <= 0) = NaN;

    % model remapped to UVP bins, 100-500 m
    phi_m     = squeeze(Y_daily(:,:,id_mod));
    phi_m_uvp = phi_m * remap_frac';
    S_m = bsxfun(@rdivide, phi_m_uvp, dw_uvp_mm) * 1e6;
    S_m(S_m <= 0) = NaN;

    % row a: UVP
    ax1 = subplot(2, n_sel, m);
    imagesc(x_log, z_uvp, log10(S_u));
    set(ax1, 'YDir','reverse','CLim',[cmin cmax],'FontSize',5,'Color','white', ...
        'XLim',xlim_log,'YLim',[0 510], ...
        'XTick',xtick_pos,'XTickLabel',{});
    if m == 1
        ylabel('Depth (m)','FontSize',6);
    else
        set(ax1,'YTickLabel',{});
    end
    title(date_lbls{m},'FontSize',6,'FontWeight','normal');

    % row b: model
    ax2 = subplot(2, n_sel, n_sel + m);
    imagesc(x_log, z_mod, log10(S_m));
    set(ax2, 'YDir','reverse','CLim',[cmin cmax],'FontSize',5,'Color','white', ...
        'XLim',xlim_log,'YLim',[0 510], ...
        'XTick',xtick_pos,'XTickLabel',{});
    if m == 1
        ylabel('Depth (m)','FontSize',6);
        set(ax2,'XTickLabel',xtick_lbl);
        xlabel('D (mm)','FontSize',6);
    else
        set(ax2,'YTickLabel',{});
    end
end

annotation('textbox',[0.005 0.90 0.04 0.05], ...
    'String','a)','EdgeColor','none','FontSize',8,'FontWeight','bold');
annotation('textbox',[0.005 0.42 0.04 0.05], ...
    'String','b)','EdgeColor','none','FontSize',8,'FontWeight','bold');
annotation('textbox',[0.04 0.90 0.08 0.05], ...
    'String','UVP 0-500 m','EdgeColor','none','FontSize',7);
annotation('textbox',[0.04 0.42 0.12 0.05], ...
    'String',sprintf('Model (\\alpha=%.3f, bc=%.3f, r_0=%.3f)', ...
    cfg.alpha, bc_scale, cfg.microbe_r0), ...
    'EdgeColor','none','FontSize',6);

cb = colorbar('Position',[0.945 0.06 0.012 0.88]);
cb.Label.String = 'Particle Volume Spectrum (ppmV mm^{-1})';
cb.Label.FontSize = 6;
set(cb,'Ticks',[-1 0 1],'FontSize',6);
cb.TickLabels = {'10^{-1}','10^0','10^1'};

saveas(gcf, fullfile(fig_dir, 'cast_spectrum_inverse_fit.png'));
fprintf('Saved cast_spectrum_inverse_fit.png\n');

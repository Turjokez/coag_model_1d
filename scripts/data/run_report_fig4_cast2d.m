% run_report_fig4_cast2d.m
%
% Cast-by-cast 2D spectrum, Siegel (2025) Fig 2a style.
%
% Row a: UVP  (0-500 m, 0.1-10 mm)  -- same dates and range as Siegel
% Row b: Model (100-500 m, 0.1-10 mm) -- model remapped to UVP bins
%                                        above 100m is blank (no model there)
%                                        above 2mm model shows zero (no large particles)
%
% Selected dates match Siegel Fig 2a: 9 panels spread across cruise.
%
% BC + run logic same as run_report_fig1_mismatch.m:
%   flux BC, Y reset to zero before final run (26-day transient).
%
% Saves: docs/figures/report_fig4_cast2d.png

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

% ---------------------------------------------------------------
% 1. Config + grid
% ---------------------------------------------------------------
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);
z_centers = col_grid.z_centers;

k_bc   = 2;
k_plot = 2:10;     % model layers 2-10 (centers 75-475 m, i.e. 100-500 m zone)
z_mod  = z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 80;
n_z           = col_grid.n_z;
dz            = col_grid.dz;

% ---------------------------------------------------------------
% 2. BC + UVP
% ---------------------------------------------------------------
cfg = cfg_base();
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;   % [n_days x n_sec]
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp = ColumnSimulation(cfg, col_grid, prof);
n_sec   = cfg.n_sections;
d_cm    = sim_tmp.size_grid.dcomb(:)';
w_bin   = (66 * d_cm .^ 0.62);

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
% Above 2000 um data can be sparse but we include it for visual comparison.
mask_uvp  = uvpd.d_um >= 100 & uvpd.d_um <= 10000;
d_uvp_mm  = uvpd.d_um(mask_uvp) / 1000;
dw_uvp_mm = uvpd.dw(mask_uvp)   / 1000;
n_uvp     = sum(mask_uvp);

% UVP depth: 0-500 m (full column, Siegel style)
mask_z_uvp = uvpd.depth_m >= 0   & uvpd.depth_m <= 510;
z_uvp      = uvpd.depth_m(mask_z_uvp);

% ---------------------------------------------------------------
% Remap: model bins -> UVP bins
%   remap_frac(j,k) = fraction of model bin k in UVP bin j
%   phi_model_on_uvp(j) = sum_k[ phi_mod(k) * remap_frac(j,k) ]
%   S_model(j) [ppmV/mm] = phi_model_on_uvp(j) / dw_uvp_mm(j) * 1e6
%
% UVP bins above 2 mm get zero from model (model has no particles that large).
% ---------------------------------------------------------------
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
% 3. Spinup (flux BC)
% ---------------------------------------------------------------
fprintf('Running spinup...\n');
sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);

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

% ---------------------------------------------------------------
% 4. Final run: flux BC, Y reset (same as fig 1)
% ---------------------------------------------------------------
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);
Y_daily = zeros(numel(k_plot), n_sec, n_days);

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
    end
    Y_daily(:,:,i_day) = (Y(k_plot,:) + Yfp(k_plot,:));
end
fprintf('Model run complete\n');

% ---------------------------------------------------------------
% 5. Select dates (Siegel Fig 2a: 9 dates across cruise)
% ---------------------------------------------------------------
% Dates 05-05, 05-07, 05-12 were transit or early-arrival casts in a
% different water mass. From 05-14 onward the ship was settled on-station.
% Use 9 on-station dates spread across the cruise.
sel_dates = [20210514, 20210515, 20210516, 20210519, ...
             20210521, 20210522, 20210525, 20210527, 20210529];

[~, ia, ib] = intersect(bc.dates, uvpd.dates);   % cast days in both
cast_dates   = bc.dates(ia);                      % model day indices
uvp_cast_idx = ib;

% find selected dates within cast days
[~, isel_c, ~] = intersect(cast_dates, sel_dates);
n_sel = numel(isel_c);
if n_sel == 0
    error('None of the selected dates found in cast days. Check date format.');
end
fprintf('Selected %d dates\n', n_sel);

% build date labels: YYYYMMDD -> 'MM-DD'
date_lbls = cell(n_sel, 1);
for m = 1:n_sel
    ds = num2str(cast_dates(isel_c(m)));   % '20210505'
    date_lbls{m} = [ds(5:6) '-' ds(7:8)];
end

% ---------------------------------------------------------------
% 6. Figure: 2 rows x n_sel columns
%    Row a: UVP  (0-500 m, 0.1-10 mm)
%    Row b: Model (100-500 m remapped to UVP bins, blank above 100 m)
% ---------------------------------------------------------------
cmin = -1;   cmax = 1;   % log10(ppmV/mm), same as Siegel
xlim_log  = [log10(0.09) log10(12)];   % 0.1-10 mm (Siegel range)
xtick_pos = log10([0.1 1.0 10]);
xtick_lbl = {'0.1', '1', '10'};

x_log = log10(d_uvp_mm);   % shared x-axis (UVP bins)

fig_w = max(n_sel * 2.0 + 1.5, 16);
figure('Units','centimeters','Position',[1 1 fig_w 10],'Color','white');
colormap(jet);

for m = 1:n_sel
    id_cast = isel_c(m);           % index into cast_dates
    id_mod  = ia(id_cast);         % index into bc.dates (model day)
    id_uvp  = uvp_cast_idx(id_cast); % index into uvpd.dates

    % --- UVP spectral density [ppmV/mm], 0-500 m ---
    phi_all = squeeze(uvpd.phi(id_uvp, :, :));   % [n_depth, n_bins]
    phi_u   = phi_all(mask_z_uvp, mask_uvp);      % [n_z_uvp, n_uvp]
    phi_u(isnan(phi_u)) = 0;
    S_u = bsxfun(@rdivide, phi_u, dw_uvp_mm) * 1e6;
    S_u(S_u <= 0) = NaN;

    % --- Model remapped to UVP bins, 100-500 m ---
    phi_m     = squeeze(Y_daily(:,:,id_mod));      % [n_z x n_sec]
    phi_m_uvp = phi_m * remap_frac';               % [n_z x n_uvp]
    S_m = bsxfun(@rdivide, phi_m_uvp, dw_uvp_mm) * 1e6;
    S_m(S_m <= 0) = NaN;

    % --- Row a: UVP (0-500 m) ---
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

    % --- Row b: Model (100-500 m, blank above) ---
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

% row labels
annotation('textbox',[0.005 0.90 0.04 0.05], ...
    'String','a)','EdgeColor','none','FontSize',8,'FontWeight','bold');
annotation('textbox',[0.005 0.42 0.04 0.05], ...
    'String','b)','EdgeColor','none','FontSize',8,'FontWeight','bold');

annotation('textbox',[0.04 0.90 0.08 0.05], ...
    'String','UVP 0-500 m','EdgeColor','none','FontSize',7);
annotation('textbox',[0.04 0.42 0.10 0.05], ...
    'String','Model 100-500 m','EdgeColor','none','FontSize',7);

% colorbar
cb = colorbar('Position',[0.945 0.06 0.012 0.88]);
cb.Label.String = 'Particle Volume Spectrum (ppmV mm^{-1})';
cb.Label.FontSize = 6;
set(cb,'Ticks',[-1 0 1],'FontSize',6);
cb.TickLabels = {'10^{-1}','10^0','10^1'};

saveas(gcf, fullfile(fig_dir, 'report_fig4_cast2d.png'));
fprintf('Saved report_fig4_cast2d.png\n');

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

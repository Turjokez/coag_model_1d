% run_report_fig1_mismatch.m
%
% Figure 1 for combined report.
% Shows where the model matches UVP and where it does not.
%
% Panel (a): Model BV vs UVP BV profile (log x-axis)
% Panel (b): Model/UVP ratio vs depth
%            Red dot = overestimate (ratio > 1)
%            Blue dot = underestimate (ratio < 1)
%            Shaded region: factor-of-2 band (0.5 to 2)
%
% Config: alpha=0.10, Da*5, zoo on, mining on, flux BC at 100 m.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir  = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- grid and turbulence ---
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
k_plot        = 2:10;          % layers 2-10 = 75-475 m
z_mod         = z_centers(k_plot);
k_bc          = 2;

% --- load UVP ---
cfg0 = cfg_base();
bc   = get_daily_bc_at_depth(uvp_file, cfg0, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp  = ColumnSimulation(cfg0, col_grid, prof);
d_cm     = sim_tmp.size_grid.dcomb(:)';
w_bin    = (66 * d_cm .^ 0.62);
mask_uvp = bc.d_model_um >= 100 & bc.d_model_um < 2000;

% UVP reference (cast-day average, 100-2000 um)
mask_uvp_raw = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z_uvp   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ~, ib]   = intersect(bc.dates, uvpd.dates);

phi_uvp = zeros(numel(k_plot), 1);
n_cast_uvp = 0;
for m = 1:numel(ib)
    phi_u = squeeze(uvpd.phi(ib(m), mask_z_uvp, mask_uvp_raw));
    if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z_uvp) - z_mod(ki)));
        phi_uvp(ki) = phi_uvp(ki) + sum(phi_u(iz,:));
    end
    n_cast_uvp = n_cast_uvp + 1;
end
phi_uvp = phi_uvp / max(n_cast_uvp, 1);
mask_ok = phi_uvp >= 0.01 * max(phi_uvp);

% --- run model: spinup ---
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

% --- final run: accumulate on cast days ---
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);
phi_mod = zeros(numel(k_plot), 1);
n_cast  = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Ytot = Y(k_plot,:) + Yfp(k_plot,:);
        phi_mod = phi_mod + sum(Ytot(:, mask_uvp), 2);
        n_cast = n_cast + 1;
    end
end
phi_mod = phi_mod / max(n_cast, 1);

% ratio
ratio = NaN(numel(k_plot), 1);
ratio(mask_ok) = phi_mod(mask_ok) ./ phi_uvp(mask_ok);

fprintf('\nDepth | Model/UVP\n');
for ki = 1:numel(k_plot)
    if mask_ok(ki)
        fprintf('  %3.0f m | %.2f\n', z_mod(ki), ratio(ki));
    end
end

% Arial font for all elements
set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

fs = 7;   % single font size for everything

% --- figure ---
figure('Units','centimeters','Position',[2 2 14 10],'Color','white');

% Panel (a): BV profile — only depths where UVP has valid data
ax1 = subplot(1,2,1);
hold on;
plot(phi_uvp(mask_ok), z_mod(mask_ok), 'b-o', 'MarkerSize',4,'LineWidth',1.2);
plot(phi_mod(mask_ok), z_mod(mask_ok), 'k-o', 'MarkerSize',4,'LineWidth',1.2);
set(ax1,'YDir','reverse','XScale','log','YLim',[60 510],'FontSize',fs,'Box','off');
xlabel('BV (m^3 m^{-3})','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);
legend({'UVP','Model'},'Location','southeast','FontSize',fs,'Box','off');
title('(a)','FontWeight','normal','FontSize',fs);

% Panel (b): ratio, color coded
ax2 = subplot(1,2,2);
hold on;

% grey band: factor of 2
fill([0.5 2 2 0.5],[60 60 510 510],[0.92 0.92 0.92],'EdgeColor','none', ...
    'HandleVisibility','off');

% ratio=1 line
plot([1 1],[60 510],'k:','LineWidth',0.8,'HandleVisibility','off');

% split into over/under and plot as two groups (clean legend)
ratio_ok = ratio;
ratio_ok(~mask_ok) = NaN;

mask_over  = ratio_ok >= 1;
mask_under = ratio_ok <  1;

% connecting line first (behind dots)
plot(ratio_ok, z_mod, 'k-','LineWidth',0.8,'HandleVisibility','off');

% dots: underestimate (all cases here are blue)
plot(ratio_ok(mask_under), z_mod(mask_under), 'o', ...
    'Color',[0.1 0.3 0.8],'MarkerFaceColor',[0.1 0.3 0.8],'MarkerSize',6, ...
    'DisplayName','underestimate');

% dots: overestimate (red, if any)
if any(mask_over)
    plot(ratio_ok(mask_over), z_mod(mask_over), 'o', ...
        'Color',[0.8 0.1 0.1],'MarkerFaceColor',[0.8 0.1 0.1],'MarkerSize',6, ...
        'DisplayName','overestimate');
end

set(ax2,'YDir','reverse','XScale','log','YLim',[60 510],...
    'XLim',[0.05 5],'FontSize',fs,'Box','off');
xlabel('Model / UVP','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);
title('(b)','FontWeight','normal','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');

saveas(gcf, fullfile(fig_dir, 'report_fig1_mismatch.png'));
fprintf('\nSaved report_fig1_mismatch.png\n');

% --- config ---
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

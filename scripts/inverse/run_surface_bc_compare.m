% run_surface_bc_compare.m
% Run the model with surface BC (k=1, z=25m) from UVP surface concentration.
% Compare to 100m BC run and UVP observations.
%
% Surface BC uses get_daily_surface_phi (depth <= 5m UVP data).
% 100m BC uses get_daily_bc_at_depth (same as before).
%
% Both runs use: alpha=0.093, r0=0.014 (best fit from 3-param inverse).
% bc_scale = 1.0 for both (no manual correction).

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));
addpath(fullfile(script_dir, '..', 'data'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);
dz       = col_grid.dz;

% best-fit config
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.r_to_rg        = 1.6;
cfg.alpha          = 0.093;
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_zoo     = true;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.enable_mining  = true;
cfg.enable_microbe = true;
cfg.microbe_r0     = 0.014;

dt            = 0.25;
steps_per_day = round(1/dt);
spinup_tol    = 0.01;
max_cycles    = 80;

% get sinking speed for BC flux
sim_tmp  = ColumnSimulation(cfg, col_grid, prof);
d_cm     = sim_tmp.size_grid.dcomb(:)';
w_bin    = 66 * d_cm .^ 0.62;   % kriest_8 [m/day]
d_um     = d_cm * 1e4;
bin_mask = d_um >= 100 & d_um < 2000;
clear sim_tmp;

% --- Run 1: surface BC (k=1, z=25m) ---
fprintf('=== Run 1: surface BC (UVP depth<=5m) ===\n');
surf = get_daily_surface_phi(uvp_file, cfg, col_grid);
phi_surf_daily = surf.phi;
n_days_surf    = surf.n_days;
fprintf('  Surface UVP: %d days, %d model bins\n', n_days_surf, cfg.n_sections);

% check surface concentrations
bv_surf_day = sum(phi_surf_daily(:, bin_mask), 2);
fprintf('  Surface BV (100-2000 um): mean=%.2e, max=%.2e cm3/cm3\n', ...
    mean(bv_surf_day), max(bv_surf_day));

Y_surf   = zeros(col_grid.n_z, cfg.n_sections);
Yfp_surf = zeros(col_grid.n_z, cfg.n_sections);
sim_surf = ColumnSimulation(cfg, col_grid, prof);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y_surf + Yfp_surf, 2));
    for i_day = 1:n_days_surf
        sim_surf.rhs.profile.eps = keps_day.eps(:, min(i_day, size(keps_day.eps,2)));
        % flux = w * phi_surface / dz  [same unit as model concentration per step]
        flux_src = dt * (w_bin .* phi_surf_daily(i_day,:)) / dz;
        for i_step = 1:steps_per_day
            Y_surf(1,:) = Y_surf(1,:) + flux_src;   % inject at k=1
            [Y_surf, Yfp_surf] = sim_surf.rhs.stepY(Y_surf, dt, Yfp_surf);
        end
    end
    phi_after = mean(sum(Y_surf + Yfp_surf, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('  Converged at cycle %d\n', icyc); break
    end
end

bv_surf_col = sum((Y_surf + Yfp_surf) .* bin_mask, 2);

% --- Run 2: 100m BC (k=2, z=75m) with bc_scale=0.42 (best fit) ---
fprintf('\n=== Run 2: 100m BC (bc_scale=0.42, best fit) ===\n');
bc100 = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 2:10);
phi_bc_100   = bc100.phi_bc_daily;
n_days_100   = bc100.n_days;

Y_100   = zeros(col_grid.n_z, cfg.n_sections);
Yfp_100 = zeros(col_grid.n_z, cfg.n_sections);
sim_100 = ColumnSimulation(cfg, col_grid, prof);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y_100 + Yfp_100, 2));
    for i_day = 1:n_days_100
        sim_100.rhs.profile.eps = keps_day.eps(:, min(i_day, size(keps_day.eps,2)));
        flux_src = dt * (w_bin .* phi_bc_100(i_day,:)) * 0.42 / dz;
        for i_step = 1:steps_per_day
            Y_100(2,:) = Y_100(2,:) + flux_src;   % inject at k=2
            [Y_100, Yfp_100] = sim_100.rhs.stepY(Y_100, dt, Yfp_100);
        end
    end
    phi_after = mean(sum(Y_100 + Yfp_100, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('  Converged at cycle %d\n', icyc); break
    end
end

bv_100_col = sum((Y_100 + Yfp_100) .* bin_mask, 2);

% --- UVP observations ---
uvp      = parse_uvp(uvp_file);
uvp_mask = uvp.d_um >= 100 & uvp.d_um < 2000;
bv_obs   = sum(uvp.phi(:, uvp_mask), 2);
dep_obs  = uvp.depth_m;

% --- ratio vs depth ---
dep_valid = dep_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);
obs_valid = bv_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);

bv_surf_interp = interp1(col_grid.z_centers, bv_surf_col, dep_valid, 'linear');
bv_100_interp  = interp1(col_grid.z_centers, bv_100_col,  dep_valid, 'linear');
ratio_surf = bv_surf_interp ./ obs_valid;
ratio_100  = bv_100_interp  ./ obs_valid;

% --- print summary ---
fprintf('\n--- Layer comparison (100-2000 um total BV) ---\n');
fprintf('  depth   UVP obs      surface BC   100m BC\n');
for k = 1:col_grid.n_z
    z = col_grid.z_centers(k);
    [~, iz] = min(abs(dep_obs - z));
    if dep_obs(iz) > z+5, continue, end
    fprintf('  %4.0f m  %.2e   %.2e   %.2e\n', z, bv_obs(iz), ...
        bv_surf_col(k), bv_100_col(k));
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 16 12],'Color','white');

subplot(1,2,1);
hold on;
plot(bv_obs, dep_obs, 'ko', 'MarkerSize',3, 'MarkerFaceColor','k', 'DisplayName','UVP obs');
plot(bv_surf_col, col_grid.z_centers, 'r-',  'LineWidth',1.2, 'DisplayName','surface BC (z=0)');
plot(bv_100_col,  col_grid.z_centers, 'b--', 'LineWidth',1.2, 'DisplayName','100m BC (bc=0.42)');
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on', ...
    'YLim',[0 1000],'XLim',[1e-7 1e-4]);
xlabel('BV 100-2000 \mum (m^3 m^{-3})','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('profile: surface vs 100m BC','FontWeight','normal','FontSize',fs);

subplot(1,2,2);
hold on;
fill([0.1 100 100 0.1],[0 0 75 75],       [0.85 0.95 0.85],'EdgeColor','none','HandleVisibility','off');
fill([0.1 100 100 0.1],[75 75 200 200],   [0.95 0.92 0.80],'EdgeColor','none','HandleVisibility','off');
fill([0.1 100 100 0.1],[200 200 1000 1000],[0.95 0.85 0.85],'EdgeColor','none','HandleVisibility','off');
plot([1 1],[0 1000],'k:','LineWidth',0.8,'HandleVisibility','off');
plot(ratio_surf, dep_valid, 'r-',  'LineWidth',1.2,'DisplayName','surface BC');
plot(ratio_100,  dep_valid, 'b--', 'LineWidth',1.2,'DisplayName','100m BC');
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on','Layer','top', ...
    'YLim',[0 1000],'XLim',[0.1 100]);
xlabel('model / obs ratio','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('ratio (1 = perfect)','FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir,'..','..','docs','figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir,'surface_bc_compare.png'));
fprintf('\nSaved surface_bc_compare.png\n');

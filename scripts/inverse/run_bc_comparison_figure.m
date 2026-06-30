% run_bc_comparison_figure.m
% Clean comparison: surface BC fit vs 100m BC fit vs UVP.
%
% Surface BC fit: alpha=0.099, bc_scale=0.057, r0=0.020 (from run_inverse_surface)
% 100m BC fit:    alpha=0.093, bc_scale=0.420, r0=0.014 (from run_inverse_3param)

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

cfg = SimulationConfig();
cfg.n_sections    = 30;
cfg.sinking_law   = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.r_to_rg       = 1.6;
cfg.enable_coag   = true;
cfg.enable_disagg = true;
cfg.disagg_mode   = 'operator_split';
cfg.disagg_dmax_A = 9.39e-6 * 5;
cfg.enable_zoo    = true;
cfg.zoo_c         = 0.025;
cfg.zoo_s         = 1.3e-5;
cfg.zoo_p         = 0.5;
cfg.zoo_ic        = 7;
cfg.enable_mining = true;
cfg.enable_microbe = false;

% bin mask 100-2000 um
sim_tmp  = ColumnSimulation(cfg, col_grid, prof);
d_um     = sim_tmp.size_grid.dcomb(:)' * 1e4;
bin_mask = d_um >= 100 & d_um < 2000;
clear sim_tmp;

all_depths = col_grid.z_centers;

% --- Run 1: surface BC fit ---
fprintf('Running surface BC fit...\n');
surf = get_daily_surface_phi(uvp_file, cfg, col_grid);
phi_surf = fwd_column_surface([0.099, 0.057, 0.020], all_depths, col_grid, ...
    keps_day, prof, surf.phi, surf.n_days, cfg);
bv_surf = sum(phi_surf .* bin_mask, 2);
fprintf('  done.\n');

% --- Run 2: 100m BC fit ---
fprintf('Running 100m BC fit...\n');
bc100 = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 2:10);
phi_100 = fwd_column_v2([0.093, 0.420, 0.014], all_depths, col_grid, ...
    keps_day, prof, bc100.phi_bc_daily, bc100.n_days, cfg);
bv_100 = sum(phi_100 .* bin_mask, 2);
fprintf('  done.\n');

% --- UVP obs ---
uvp      = parse_uvp(uvp_file);
uvp_mask = uvp.d_um >= 100 & uvp.d_um < 2000;
bv_obs   = sum(uvp.phi(:, uvp_mask), 2);
dep_obs  = uvp.depth_m;

dep_v = dep_obs(dep_obs>=25 & dep_obs<=975 & bv_obs>0);
obs_v = bv_obs(dep_obs>=25 & dep_obs<=975 & bv_obs>0);
r_surf = interp1(all_depths, bv_surf, dep_v, 'linear') ./ obs_v;
r_100  = interp1(all_depths, bv_100,  dep_v, 'linear') ./ obs_v;

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 16 12],'Color','white');

subplot(1,2,1);
hold on;
plot(bv_obs,  dep_obs,    'ko', 'MarkerSize',3,'MarkerFaceColor','k','DisplayName','UVP obs');
plot(bv_surf, all_depths, 'r-', 'LineWidth',1.4, ...
    'DisplayName','surface BC (\alpha=0.099, bc=0.057, r_0=0.020)');
plot(bv_100,  all_depths, 'b--','LineWidth',1.4, ...
    'DisplayName','100m BC (\alpha=0.093, bc=0.420, r_0=0.014)');
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on', ...
    'YLim',[0 1000],'XLim',[1e-7 1e-4]);
xlabel('BV 100-2000 \mum (m^3 m^{-3})','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('profile','FontWeight','normal','FontSize',fs);

subplot(1,2,2);
hold on;
fill([0.1 10 10 0.1],[0 0 1000 1000],[0.95 0.85 0.85],'EdgeColor','none','HandleVisibility','off');
fill([0.5 2 2 0.5],  [0 0 1000 1000],[0.85 0.95 0.85],'EdgeColor','none','HandleVisibility','off');
plot([1 1],[0 1000],'k:','LineWidth',0.8,'HandleVisibility','off');
plot(r_surf, dep_v, 'r-', 'LineWidth',1.4,'DisplayName','surface BC fit');
plot(r_100,  dep_v, 'b--','LineWidth',1.4,'DisplayName','100m BC fit');
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on','Layer','top', ...
    'YLim',[0 1000],'XLim',[0.1 10]);
xlabel('model / obs ratio','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('ratio (1 = perfect)','FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir,'..','..','docs','figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir,'bc_comparison_final.png'));
fprintf('Saved bc_comparison_final.png\n');

% print summary table
fprintf('\n--- Summary ---\n');
fprintf('Method        alpha   bc_scale   r0       J\n');
fprintf('Surface BC    0.099   0.057      0.020    1.579\n');
fprintf('100m BC       0.093   0.420      0.014    0.241\n');
fprintf('\nNote: bc_scale interpretation\n');
fprintf('  Surface (5m):  6%% of UVP surface flux enters as sinking flux\n');
fprintf('  100m:         42%% of UVP 100m conc enters as sinking flux\n');

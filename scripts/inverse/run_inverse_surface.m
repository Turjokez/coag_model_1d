% run_inverse_surface.m
% Fit [alpha, bc_scale, r0] to EXPORTS-NA UVP observations.
% Surface BC: top 5m UVP data with power-law fill, injected at k=1.
%
% alpha    : stickiness
% bc_scale : scale on surface flux (1 = raw UVP flux, <1 = reduce)
% r0       : microbial remineralization rate [day^-1]
%
% Obs: total BV at 125, 325, 475 m (100-2000 um)
% Method: fminsearch, log-space cost + prior penalty
%
% Expected runtime: ~15-20 min

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

cfg_base = SimulationConfig();
cfg_base.n_sections    = 30;
cfg_base.sinking_law   = 'kriest_8';
cfg_base.ds_kernel_mode = 'sinking_law';
cfg_base.r_to_rg       = 1.6;
cfg_base.alpha         = 0.10;
cfg_base.enable_coag   = true;
cfg_base.enable_disagg = true;
cfg_base.disagg_mode   = 'operator_split';
cfg_base.disagg_dmax_A = 9.39e-6 * 5;
cfg_base.enable_zoo    = true;
cfg_base.zoo_c         = 0.025;
cfg_base.zoo_s         = 1.3e-5;
cfg_base.zoo_p         = 0.5;
cfg_base.zoo_ic        = 7;
cfg_base.enable_microbe = false;
cfg_base.enable_mining  = true;

% surface BC from UVP top 5m (with power-law fill)
surf = get_daily_surface_phi(uvp_file, cfg_base, col_grid);
phi_surf_daily = surf.phi;
n_days         = surf.n_days;
fprintf('Surface BC: %d days, surface BV mean=%.2e\n', n_days, ...
    mean(sum(phi_surf_daily(:, 10:end), 2)));

% obs depths
obs_depths = [125, 325, 475];
obs = load_uvp_obs(uvp_file, obs_depths);
fprintf('UVP obs: %d depths\n', numel(obs_depths));

% priors — bc_scale center at 0.15 (surface flux is ~6x 100m flux)
prior.alpha    = 0.10;  sig_alpha = 1.5;
prior.bc_scale = 0.15;  sig_bc    = 1.0;   % expect large reduction from surface
prior.r0       = 0.01;  sig_r0    = 1.0;

sigma_log = 0.5;

cost = @(p) cost_surface(p, obs, obs_depths, col_grid, keps_day, prof, ...
    phi_surf_daily, n_days, cfg_base, prior, sigma_log, sig_alpha, sig_bc, sig_r0);

% starting point
p0 = [0.05, 0.10, 0.02];
fprintf('\n--- Starting surface inverse fit ---\n');
fprintf('p0: alpha=%.3f  bc_scale=%.3f  r0=%.4f\n\n', p0(1), p0(2), p0(3));

opts = optimset('TolX',1e-3,'TolFun',0.1,'MaxFunEvals',300);
tic;
[p_fit, J_fit] = fminsearch(cost, p0, opts);
t_total = toc;

fprintf('\n--- Done in %.1f min ---\n', t_total/60);
fprintf('Best fit: alpha=%.4f  bc_scale=%.4f  r0=%.5f  J=%.3f\n', ...
    p_fit(1), p_fit(2), p_fit(3), J_fit);

% evaluate
phi_fit = fwd_column_surface(p_fit, obs_depths, col_grid, keps_day, prof, ...
    phi_surf_daily, n_days, cfg_base);
bv_fit  = sum(phi_fit, 2);

phi_p0 = fwd_column_surface(p0, obs_depths, col_grid, keps_day, prof, ...
    phi_surf_daily, n_days, cfg_base);
bv_p0  = sum(phi_p0, 2);

fprintf('\nComparison:\n  depth   obs          start        fit\n');
for id = 1:numel(obs_depths)
    fprintf('  %3d m   %.3e   %.3e   %.3e\n', ...
        obs_depths(id), obs.bv_total(id), bv_p0(id), bv_fit(id));
end

% --- figure: fit at 3 obs depths ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 9],'Color','white');
hold on;
plot(obs.bv_total, obs_depths, 'ko', 'MarkerSize',5, 'MarkerFaceColor','k', ...
    'DisplayName','UVP obs');
plot(bv_p0,  obs_depths, 'b--', 'LineWidth',1.2, ...
    'DisplayName', sprintf('start (a=%.2f, bc=%.2f, r0=%.3f)', p0(1), p0(2), p0(3)));
plot(bv_fit, obs_depths, 'r-',  'LineWidth',1.2, ...
    'DisplayName', sprintf('fit (a=%.3f, bc=%.3f, r0=%.4f)', p_fit(1), p_fit(2), p_fit(3)));
set(gca,'YDir','reverse','FontSize',fs,'Box','on','XScale','log', ...
    'YLim',[0 600],'XLim',[5e-7 5e-5]);
xlabel('total BV 100-2000 \mum (m^3 m^{-3})','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('surface BC inverse fit','FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir,'..','..','docs','figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir,'inverse_fit_surface.png'));
fprintf('\nSaved inverse_fit_surface.png\n');

% --- full depth comparison ---
fprintf('\nRunning full depth comparison...\n');

% run full column for fit and prior (alpha=0.10, bc=1 reference)
all_depths = col_grid.z_centers;

% fit
phi_full_fit = fwd_column_surface(p_fit, all_depths, col_grid, keps_day, prof, ...
    phi_surf_daily, n_days, cfg_base);

% reference: alpha=0.10, bc_scale=1.0, r0=0 (no power-law fix comparison)
p_ref = [0.10, 1.00, 0.00];
phi_full_ref = fwd_column_surface(p_ref, all_depths, col_grid, keps_day, prof, ...
    phi_surf_daily, n_days, cfg_base);

% also load 100m best fit for comparison
bc100 = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, 2:10);
phi_bc100 = fwd_column_v2([0.093, 0.42, 0.014], all_depths, col_grid, keps_day, prof, ...
    bc100.phi_bc_daily, bc100.n_days, cfg_base);

sim_tmp  = ColumnSimulation(cfg_base, col_grid, prof);
d_um     = sim_tmp.size_grid.dcomb(:)' * 1e4;
bin_mask = d_um >= 100 & d_um < 2000;
clear sim_tmp;

bv_ref  = sum(phi_full_ref  .* bin_mask, 2);
bv_ffit = sum(phi_full_fit  .* bin_mask, 2);
bv_100  = sum(phi_bc100     .* bin_mask, 2);

uvp      = parse_uvp(uvp_file);
uvp_mask = uvp.d_um >= 100 & uvp.d_um < 2000;
bv_obs   = sum(uvp.phi(:, uvp_mask), 2);
dep_obs  = uvp.depth_m;

dep_v = dep_obs(dep_obs>=25 & dep_obs<=975 & bv_obs>0);
obs_v = bv_obs(dep_obs>=25 & dep_obs<=975 & bv_obs>0);

r_ref  = interp1(all_depths, bv_ref,  dep_v, 'linear') ./ obs_v;
r_ffit = interp1(all_depths, bv_ffit, dep_v, 'linear') ./ obs_v;
r_100  = interp1(all_depths, bv_100,  dep_v, 'linear') ./ obs_v;

figure('Units','centimeters','Position',[2 2 18 12],'Color','white');

subplot(1,2,1);
hold on;
plot(bv_obs,  dep_obs,    'ko','MarkerSize',3,'MarkerFaceColor','k','DisplayName','UVP obs');
plot(bv_ref,  all_depths, 'b--','LineWidth',1.2,'DisplayName','surf BC (bc=1, prior)');
plot(bv_ffit, all_depths, 'r-', 'LineWidth',1.2,'DisplayName', ...
    sprintf('surf BC fit (bc=%.2f)', p_fit(2)));
plot(bv_100,  all_depths, 'g:', 'LineWidth',1.2,'DisplayName','100m BC fit');
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on', ...
    'YLim',[0 1000],'XLim',[1e-7 1e-4]);
xlabel('BV 100-2000 \mum (m^3 m^{-3})','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('profile','FontWeight','normal','FontSize',fs);

subplot(1,2,2);
hold on;
fill([0.1 100 100 0.1],[0 0 75 75],       [0.85 0.95 0.85],'EdgeColor','none','HandleVisibility','off');
fill([0.1 100 100 0.1],[75 75 200 200],   [0.95 0.92 0.80],'EdgeColor','none','HandleVisibility','off');
fill([0.1 100 100 0.1],[200 200 1000 1000],[0.95 0.85 0.85],'EdgeColor','none','HandleVisibility','off');
plot([1 1],[0 1000],'k:','LineWidth',0.8,'HandleVisibility','off');
plot(r_ref,  dep_v, 'b--','LineWidth',1.2,'DisplayName','surf BC (bc=1)');
plot(r_ffit, dep_v, 'r-', 'LineWidth',1.2,'DisplayName','surf BC fit');
plot(r_100,  dep_v, 'g:', 'LineWidth',1.2,'DisplayName','100m BC fit');
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on','Layer','top', ...
    'YLim',[0 1000],'XLim',[0.1 100]);
xlabel('model / obs ratio','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('ratio (1 = perfect)','FontWeight','normal','FontSize',fs);

saveas(gcf, fullfile(fig_dir,'inverse_surface_full_depth.png'));
fprintf('Saved inverse_surface_full_depth.png\n');


% ---- local cost function ----
function J = cost_surface(p, obs, obs_depths, col_grid, keps_day, prof, ...
    phi_surf_daily, n_days, cfg_base, prior, sigma_log, sig_alpha, sig_bc, sig_r0)

alpha_try = p(1); bc_sc = p(2); r0 = p(3);

if alpha_try<=0 || alpha_try>2 || bc_sc<=0 || bc_sc>2 || r0<0 || r0>1
    J = 1e6; return
end

phi_mod = fwd_column_surface(p, obs_depths, col_grid, keps_day, prof, ...
    phi_surf_daily, n_days, cfg_base);
bv_mod  = sum(phi_mod, 2);

J_data = 0;
for id = 1:numel(obs_depths)
    if bv_mod(id) <= 0, J = 1e6; return, end
    J_data = J_data + ((log(bv_mod(id)) - log(obs.bv_total(id))) / sigma_log)^2;
end

J_prior = ((log(alpha_try)   - log(prior.alpha))    / sig_alpha)^2 + ...
          ((log(bc_sc)        - log(prior.bc_scale)) / sig_bc)^2    + ...
          ((log(r0+1e-8)      - log(prior.r0))       / sig_r0)^2;

J = J_data + J_prior;
fprintf('  a=%.3f bc=%.3f r0=%.4f  Jd=%.2f Jp=%.2f J=%.2f\n', ...
    alpha_try, bc_sc, r0, J_data, J_prior, J);
end

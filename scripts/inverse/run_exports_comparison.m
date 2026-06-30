% run_exports_comparison.m
% Compare best-fit model vs EXPORTS-NA UVP observations.
%
% Best fit (June 25 2026): alpha=0.093, bc_scale=0.42, r0=0.014
% Prior (old default):     alpha=0.10,  bc_scale=1.00, r0=0.000
%
% Figure 1: BV profile vs UVP obs + model/obs ratio
% Figure 2: BV flux profile (Th-234 scaled at 100 m)

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
z_c      = col_grid.z_centers;

cfg_base = SimulationConfig();
cfg_base.n_sections     = 30;
cfg_base.sinking_law    = 'kriest_8';
cfg_base.ds_kernel_mode = 'sinking_law';
cfg_base.r_to_rg        = 1.6;
cfg_base.enable_coag    = true;
cfg_base.enable_disagg  = true;
cfg_base.disagg_mode    = 'operator_split';
cfg_base.disagg_dmax_A  = 9.39e-6 * 5;
cfg_base.enable_zoo     = true;
cfg_base.zoo_c          = 0.025;
cfg_base.zoo_s          = 1.3e-5;
cfg_base.zoo_p          = 0.5;
cfg_base.zoo_ic         = 7;
cfg_base.enable_mining  = true;

bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, 2:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

% model bin range 100-2000 um (UVP range)
sim_tmp  = ColumnSimulation(cfg_base, col_grid, prof);
d_um     = sim_tmp.size_grid.dcomb(:)' * 1e4;
bin_mask = d_um >= 100 & d_um < 2000;
clear sim_tmp;

% --- two runs ---
% [alpha, bc_scale, r0, enable_microbe]
runs = {[0.10,  1.00, 0.000, false], ...   % prior
        [0.093, 0.42, 0.014, true  ]};     % best fit
labels = {'prior (\alpha=0.10, bc=1.0, r_0=0)', ...
          'fit  (\alpha=0.093, bc=0.42, r_0=0.014)'};
colors = {'b', 'r'};
lstyle = {'--', '-'};

dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;
spinup_tol    = 0.01;
max_cycles    = 80;

BV_store  = zeros(col_grid.n_z, 2);   % total BV in 100-2000 um
FLX_store = zeros(col_grid.n_z, 2);   % BV flux

for ia = 1:2
    p   = runs{ia};
    cfg = copy(cfg_base);
    cfg.alpha          = p(1);
    cfg.enable_microbe = p(4);
    cfg.microbe_r0     = p(3);

    sim   = ColumnSimulation(cfg, col_grid, prof);
    w_bin = 66 * sim.size_grid.dcomb(:)' .^ 0.62;
    bc_sc = p(2);

    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);

    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) * bc_sc / dz;
            for i_step = 1:steps_per_day
                Y(k_bc,:) = Y(k_bc,:) + flux_src;
                [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
            end
        end
        phi_after = mean(sum(Y + Yfp, 2));
        if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
            fprintf('%s: converged cycle %d\n', labels{ia}, icyc);
            break
        end
    end

    BV_store(:, ia)  = sum((Y + Yfp) .* bin_mask, 2);
    FLX_store(:, ia) = sum((Y + Yfp) .* w_bin,    2);
end

% --- UVP observations ---
uvp      = parse_uvp(uvp_file);
uvp_mask = uvp.d_um >= 100 & uvp.d_um < 2000;
bv_obs   = sum(uvp.phi(:, uvp_mask), 2);
dep_obs  = uvp.depth_m;

dep_valid = dep_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);
obs_valid = bv_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);

% --- Figure 1: BV profile + ratio ---
fs = 7;
figure('Units','centimeters','Position',[2 2 16 12],'Color','white');

subplot(1,2,1);
hold on;
plot(bv_obs, dep_obs, 'ko', 'MarkerSize',3, 'MarkerFaceColor','k', 'DisplayName','UVP obs');
for ia = 1:2
    plot(BV_store(:,ia), z_c, [colors{ia} lstyle{ia}], ...
        'LineWidth',1.2, 'DisplayName', labels{ia});
end
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on', ...
    'YLim',[0 1000],'XLim',[1e-7 1e-4]);
xlabel('BV 100-2000 \mum (m^3 m^{-3})','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('BV profile vs UVP','FontWeight','normal','FontSize',fs);

subplot(1,2,2);
hold on;
for ia = 1:2
    mod_interp = interp1(z_c, BV_store(:,ia), dep_valid, 'linear');
    ratio = mod_interp ./ obs_valid;
    plot(ratio, dep_valid, [colors{ia} lstyle{ia}], ...
        'LineWidth',1.2, 'DisplayName', labels{ia});
end
xline(1,'k:','LineWidth',0.8,'HandleVisibility','off');
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on', ...
    'YLim',[0 1000],'XLim',[0.1 100]);
xlabel('model / obs ratio','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('ratio (1 = perfect)','FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir,'..','..','docs','figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
exportgraphics(gcf, fullfile(fig_dir,'exports_bv_comparison.png'), 'Resolution',200);
fprintf('\nSaved exports_bv_comparison.png\n');

% --- print summary ---
fprintf('\n--- Model / UVP ratio at key depths (fit config) ---\n');
for zz = [125, 225, 325, 475, 725]
    [~, ku] = min(abs(dep_valid - zz));
    [~, km] = min(abs(z_c - zz));
    r = BV_store(km, 2) / obs_valid(ku);
    fprintf('  z = %4.0f m: model/obs = %.2f\n', z_c(km), r);
end

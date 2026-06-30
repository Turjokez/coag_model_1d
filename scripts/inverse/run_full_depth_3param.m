% run_full_depth_3param.m
% Full depth profile comparison: prior vs 3-param best fit vs UVP obs.
% Best fit from run_inverse_3param: alpha=0.093, bc_scale=0.42, r0=0.014

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

cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.r_to_rg        = 1.6;
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

k_plot = 2:10;
bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;
spinup_tol    = 0.01;
max_cycles    = 80;

% model bins in 100-2000 um range
sim_tmp  = ColumnSimulation(cfg, col_grid, prof);
d_um     = sim_tmp.size_grid.dcomb(:)' * 1e4;
bin_mask = d_um >= 100 & d_um < 2000;
clear sim_tmp;

% configs: [alpha, bc_scale, r0, enable_microbe]
runs = {[0.10, 1.00, 0.000, false], ...   % prior (no microbe, full BC)
        [0.093, 0.42, 0.014, true]};      % 3-param best fit

labels = {'prior (\alpha=0.10, bc=1.0, r0=0)', ...
          'fit (\alpha=0.093, bc=0.42, r0=0.014)'};
colors = {'b', 'r'};
lines  = {'--', '-'};

Y_store = cell(2,1);

for ia = 1:2
    p = runs{ia};
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
            fprintf('%s: converged cycle %d\n', labels{ia}, icyc); break
        end
    end
    Y_store{ia} = sum((Y + Yfp) .* bin_mask, 2);
end

% UVP obs at all depths
uvp     = parse_uvp(uvp_file);
uvp_mask = uvp.d_um >= 100 & uvp.d_um < 2000;
bv_obs  = sum(uvp.phi(:, uvp_mask), 2);
dep_obs = uvp.depth_m;

% ratio vs depth
dep_valid = dep_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);
obs_valid = bv_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);
ratio = zeros(numel(dep_valid), 2);
for ia = 1:2
    ratio(:,ia) = interp1(col_grid.z_centers, Y_store{ia}, dep_valid, 'linear') ./ obs_valid;
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 16 12],'Color','white');

subplot(1,2,1);
hold on;
plot(bv_obs, dep_obs, 'ko', 'MarkerSize',3, 'MarkerFaceColor','k', 'DisplayName','UVP obs');
for ia = 1:2
    plot(Y_store{ia}, col_grid.z_centers, [colors{ia} lines{ia}], ...
        'LineWidth',1.2, 'DisplayName', labels{ia});
end
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on', ...
    'YLim',[0 1000],'XLim',[1e-7 1e-4]);
xlabel('BV 100-2000 \mum (m^3 m^{-3})','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('profile','FontWeight','normal','FontSize',fs);

subplot(1,2,2);
hold on;
fill([0.1 100 100 0.1],[0 0 75 75],   [0.85 0.95 0.85],'EdgeColor','none','HandleVisibility','off');
fill([0.1 100 100 0.1],[75 75 200 200],[0.95 0.92 0.80],'EdgeColor','none','HandleVisibility','off');
fill([0.1 100 100 0.1],[200 200 1000 1000],[0.95 0.85 0.85],'EdgeColor','none','HandleVisibility','off');
plot([1 1],[0 1000],'k:','LineWidth',0.8,'HandleVisibility','off');
for ia = 1:2
    plot(ratio(:,ia), dep_valid, [colors{ia} lines{ia}], ...
        'LineWidth',1.2,'DisplayName',labels{ia});
end
text(0.12, 38,  'near match',       'FontSize',fs-1,'Color',[0.2 0.5 0.2]);
text(0.12, 138, 'onset',            'FontSize',fs-1,'Color',[0.6 0.4 0.0]);
text(0.12, 600, 'structural drift', 'FontSize',fs-1,'Color',[0.6 0.1 0.1]);
set(gca,'YDir','reverse','XScale','log','FontSize',fs,'Box','on','Layer','top', ...
    'YLim',[0 1000],'XLim',[0.1 100]);
xlabel('model / obs ratio','FontSize',fs);
ylabel('depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('ratio (1 = perfect)','FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir,'..','..','docs','figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir,'full_depth_3param.png'));
fprintf('Saved full_depth_3param.png\n');

% run_report_fig4_sensitivity.m
%
% Figure 4 for combined report.
% Model/UVP ratio vs depth for 4 configs on the same axes.
%
% Configs:
%   1. Best     : alpha=0.10, Da*5, zoo on,  BC=1.0
%   2. BC fix   : alpha=0.10, Da*5, zoo on,  BC*0.2
%   3. Disagg*1 : alpha=0.10, Da*1, zoo on,  BC=1.0
%   4. Zoo off  : alpha=0.10, Da*5, zoo off, BC=1.0
%
% This shows that:
%   - BC*0.2 fixes the near-surface mismatch but not the deep slope
%   - Reducing disagg or removing zoo does not improve deep ratio
%   - The deep attenuation problem (b too steep) is robust to these knobs
%
% Saves: report_fig4_sensitivity.png

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
z_mod         = z_centers(k_plot);
k_bc          = 2;

% --- load UVP reference ---
bc   = get_daily_bc_at_depth(uvp_file, cfg_make(1.0,5,true), col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp  = ColumnSimulation(cfg_make(1.0,5,true), col_grid, prof);
d_cm     = sim_tmp.size_grid.dcomb(:)';
w_bin    = (66 * d_cm .^ 0.62);
mask_uvp = bc.d_model_um >= 100 & bc.d_model_um < 2000;

mask_uvp_raw = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z_uvp   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ~, ib]   = intersect(bc.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ib)
    phi_u = squeeze(uvpd.phi(ib(m), mask_z_uvp, mask_uvp_raw));
    if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z_uvp) - z_mod(ki)));
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz,:));
    end
end
phi_uvp_ref = phi_uvp_ref / max(numel(ib), 1);
mask_ok     = phi_uvp_ref >= 0.01 * max(phi_uvp_ref);

% --- configs ---
cfglabels  = {'Best (bc=1)',  'BC×0.2',  'Disagg×1',  'Zoo off'};
bc_scales  = [1.0,             0.2,        1.0,          1.0  ];
disagg_A   = [5,               5,          1,            5    ];
zoo_flags  = [true,            true,       true,         false];
linestyles = {'k-',           'b-',       'k--',        'k:' };

n_cfg     = numel(cfglabels);
ratio_out = NaN(numel(k_plot), n_cfg);

for ic = 1:n_cfg
    fprintf('\n--- %s ---\n', cfglabels{ic});
    cfg  = cfg_make(bc_scales(ic), disagg_A(ic), zoo_flags(ic));
    sim  = ColumnSimulation(cfg, col_grid, prof);
    Y    = zeros(n_z, cfg.n_sections);
    Yfp  = zeros(n_z, cfg.n_sections);
    phi_bc_scaled = phi_bc_daily * bc_scales(ic);

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            flux_src = dt * (w_bin .* phi_bc_scaled(i_day,:)) / dz;
            for i_step = 1:steps_per_day
                Y(k_bc,:) = Y(k_bc,:) + flux_src;
                [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
            end
        end
        phi_after = mean(sum(Y + Yfp, 2));
        if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
            fprintf('  Converged at cycle %d\n', icyc); break;
        end
    end

    % collect on cast days
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    phi_mod = zeros(numel(k_plot), 1);
    n_cast  = 0;
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_scaled(i_day,:)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc,:) = Y(k_bc,:) + flux_src;
            [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            Ytot    = Y(k_plot,:) + Yfp(k_plot,:);
            phi_mod = phi_mod + sum(Ytot(:, mask_uvp), 2);
            n_cast  = n_cast + 1;
        end
    end
    phi_mod = phi_mod / max(n_cast, 1);
    ratio_out(mask_ok, ic) = phi_mod(mask_ok) ./ phi_uvp_ref(mask_ok);
end

% print table
fprintf('\nDepth  ');
for ic = 1:n_cfg, fprintf('  %-12s', cfglabels{ic}); end; fprintf('\n');
for ki = 1:numel(k_plot)
    if mask_ok(ki)
        fprintf('%5.0f m', z_mod(ki));
        for ic = 1:n_cfg, fprintf('  %12.3f', ratio_out(ki,ic)); end
        fprintf('\n');
    end
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 10],'Color','white');
hold on;

fill([0.5 2 2 0.5],[min(z_mod(mask_ok)) min(z_mod(mask_ok)) ...
    max(z_mod(mask_ok)) max(z_mod(mask_ok))], ...
    [0.92 0.92 0.92],'EdgeColor','none','HandleVisibility','off');

plot([1 1],[min(z_mod(mask_ok)) max(z_mod(mask_ok))],'k:','LineWidth',0.8, ...
    'HandleVisibility','off');

for ic = 1:n_cfg
    plot(ratio_out(:,ic), z_mod, linestyles{ic}, ...
        'MarkerSize',3,'LineWidth',1.2,'DisplayName',cfglabels{ic});
end

set(gca,'YDir','reverse','XScale','log','YLim',[60 510],'XLim',[0.05 5], ...
    'FontSize',fs,'Box','off');
xlabel('Model / UVP','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('Sensitivity','FontWeight','normal','FontSize',fs);

saveas(gcf, fullfile(fig_dir, 'report_fig4_sensitivity.png'));
fprintf('Saved report_fig4_sensitivity.png\n');

% --- config factory ---
function cfg = cfg_make(bc_scale, disagg_fac, enable_zoo)
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.disagg_dmax_A  = 9.39e-6 * disagg_fac;
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.enable_zoo     = enable_zoo;
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
% bc_scale applied outside: phi_bc_daily * bc_scale
end

% run_report_fig4_martin.m
%
% Figure 4 for combined report.
% Martin curve: model BV flux vs depth, compared to canonical power laws.
%
% F(z) = F_ref * (z / z_ref)^(-b)
%
% We plot three curves, all normalized to F(100 m):
%   1. Model (best config, flux BC)
%   2. Canonical: b = 0.858  (Martin et al. 1987)
%   3. Model fit: b = 1.72   (what our model actually gives)
%
% This is the clearest view of Problem 2:
% the model attenuates particle flux at roughly twice the observed rate.
%
% Saves: report_fig4_martin.png

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

cfg = cfg_base();
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp = ColumnSimulation(cfg, col_grid, prof);
d_cm    = sim_tmp.size_grid.dcomb(:)';
w_bin   = (66 * d_cm .^ 0.62);

% --- run model ---
fprintf('Running model...\n');
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

% BV flux at each layer, averaged over cast days
[~, ~, ib] = intersect(bc.dates, uvpd.dates);
F_bv_sum   = zeros(n_z, 1);
n_cast      = 0;
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        for k = 1:n_z
            F_bv_sum(k) = F_bv_sum(k) + sum(w_bin .* Y(k,:));
        end
        n_cast = n_cast + 1;
    end
end
F_bv = F_bv_sum / max(n_cast, 1);

% --- reference depth and normalization ---
z_ref = 100;
[~, k_ref] = min(abs(z_centers - z_ref));
[~, k500]  = min(abs(z_centers - 500));

F_norm = F_bv / F_bv(k_ref);   % normalize to 1 at z_ref

% fit b from model (100 to 500 m)
b_mod = -log(F_bv(k500) / F_bv(k_ref)) / log(z_centers(k500) / z_ref);
fprintf('Model Martin b = %.3f\n', b_mod);

% --- power law curves ---
z_line  = (50:10:600)';
b_canon = 0.858;
F_canon = (z_line / z_ref) .^ (-b_canon);   % canonical
F_fit   = (z_line / z_ref) .^ (-b_mod);     % model fit

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 10],'Color','white');
hold on;

% shaded region between the two power laws
fill([F_canon; flipud(F_fit)], [z_line; flipud(z_line)], ...
    [0.90 0.90 0.90], 'EdgeColor','none','HandleVisibility','off');

% model curve
plot(F_norm(k_ref:k500), z_centers(k_ref:k500), 'k-o', ...
    'MarkerSize',4,'LineWidth',1.4,'DisplayName','Model');

% canonical b
plot(F_canon, z_line, 'b--', 'LineWidth',1.2, ...
    'DisplayName',sprintf('b = %.3f (Martin 1987)', b_canon));

% model-fit b
plot(F_fit, z_line, 'k:', 'LineWidth',1.2, ...
    'DisplayName',sprintf('b = %.2f (model fit)', b_mod));

set(gca,'YDir','reverse','XScale','log','YLim',[80 550],'XLim',[0.01 1.2], ...
    'FontSize',fs,'Box','off');
xlabel('F(z) / F(100 m)','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('Flux attenuation','FontWeight','normal','FontSize',fs);

saveas(gcf, fullfile(fig_dir, 'report_fig4_martin.png'));
fprintf('Saved report_fig4_martin.png\n');

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

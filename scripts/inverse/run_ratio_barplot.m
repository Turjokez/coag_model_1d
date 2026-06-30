% run_ratio_barplot.m
% Model/UVP BV ratio at key depths: prior vs best fit.
%
% Simple bar chart. Dashed line at ratio=1 (perfect match).
% Depths: 125, 225, 325, 475, 725 m.
%
% Saves: docs/figs/ratio_barplot.png

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
z_c      = col_grid.z_centers;
dz       = col_grid.dz;

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

sim_tmp  = ColumnSimulation(cfg_base, col_grid, prof);
d_um     = sim_tmp.size_grid.dcomb(:)' * 1e4;
bin_mask = d_um >= 100 & d_um < 2000;
clear sim_tmp;

% two configs: [alpha, bc_scale, r0, enable_microbe]
runs   = {[0.10,  1.00, 0.000, false], ...
          [0.093, 0.42, 0.014, true ]};
rlabels = {'prior', 'best fit'};
colors  = {[0.4 0.6 0.9], [0.9 0.3 0.2]};

dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;
spinup_tol    = 0.01;
max_cycles    = 80;

BV_store = zeros(col_grid.n_z, 2);

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
            fprintf('%s: converged cycle %d\n', rlabels{ia}, icyc);
            break
        end
    end
    BV_store(:, ia) = sum((Y + Yfp) .* bin_mask, 2);
end

% --- UVP obs ---
uvp      = parse_uvp(uvp_file);
uvp_mask = uvp.d_um >= 100 & uvp.d_um < 2000;
bv_obs   = sum(uvp.phi(:, uvp_mask), 2);
dep_obs  = uvp.depth_m;
dep_valid = dep_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);
obs_valid = bv_obs(dep_obs >= 25 & dep_obs <= 975 & bv_obs > 0);

% --- compute ratios at key depths ---
plot_depths = [125, 225, 325, 475, 725];
n_dep = numel(plot_depths);
ratios = zeros(n_dep, 2);

for id = 1:n_dep
    [~, ku] = min(abs(dep_valid - plot_depths(id)));
    [~, km] = min(abs(z_c - plot_depths(id)));
    for ia = 1:2
        ratios(id, ia) = BV_store(km, ia) / obs_valid(ku);
    end
end

% print table
fprintf('\n--- Model/UVP ratio ---\n');
fprintf('  %6s  %8s  %8s\n', 'z(m)', 'prior', 'best fit');
for id = 1:n_dep
    fprintf('  %6.0f  %8.2f  %8.2f\n', plot_depths(id), ratios(id,1), ratios(id,2));
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 10],'Color','white');
hold on;
xline(1, 'k:', 'LineWidth', 0.8, 'HandleVisibility','off');
plot(ratios(:,1), plot_depths, 'o--', 'Color', [0.3 0.5 0.85], ...
    'MarkerFaceColor',[0.3 0.5 0.85], 'MarkerSize',5, 'LineWidth',1.0, ...
    'DisplayName','prior');
plot(ratios(:,2), plot_depths, 'o-', 'Color', [0.85 0.2 0.15], ...
    'MarkerFaceColor',[0.85 0.2 0.15], 'MarkerSize',5, 'LineWidth',1.2, ...
    'DisplayName','best fit');
set(gca, 'YDir','reverse', 'XScale','log', 'FontSize',fs, 'Box','on', ...
    'YLim',[100 800], 'XLim',[0.1 20], 'YTick', plot_depths);
xlabel('model / UVP BV', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','northeast', 'FontSize',fs, 'Box','off');

fig_dir = fullfile(script_dir,'..','..','docs','figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
exportgraphics(gcf, fullfile(fig_dir,'ratio_barplot.png'), 'Resolution',200);
fprintf('\nSaved ratio_barplot.png\n');

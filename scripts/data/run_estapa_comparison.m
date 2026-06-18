% run_estapa_comparison.m
%
% Compare model BV flux to Estapa STT+NBST POC flux.
%
% Model: F_bv(z) = sum_d w(d) * Y(z,d)  [m3 BV m-2 d-1]
% Trap:  flux_POC                         [mg C m-2 d-1]
%
% We compute the implied C:BV ratio = trap_POC / F_bv
% and check whether it is physically reasonable (~1e5 to 1e7 mg C m-3).
%
% Comparison depths: 75, 125, 175, 330, 500 m (STT coverage).
% Uses flux BC, best config.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path  = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file  = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
stt_file  = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'sediment_trap_estapa', 'raw', ...
    '54049d4152_EXPORTS-EXPORTSNA_JC214_STT_fluxes.sb');
nbst_file = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'sediment_trap_estapa', 'raw', ...
    '54049d4152_EXPORTS-EXPORTSNA_JC214_NBST_fluxes.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Load Estapa trap data
% ---------------------------------------------------------------
trap = load_estapa_flux(stt_file, nbst_file);
fprintf('Trap depths: '); fprintf('%d ', trap.depths); fprintf('m\n');

% comparison depths: 75, 125, 175, 330, 500 m
target_z = [75, 125, 175, 330, 500];
id_trap  = zeros(1, numel(target_z));
for i = 1:numel(target_z)
    [~, id_trap(i)] = min(abs(trap.depths - target_z(i)));
end
z_comp    = trap.depths(id_trap);
poc_comp  = trap.flux_POC(id_trap);       % mg C m-2 d-1
poc_sd    = trap.flux_POC_sd(id_trap);
src_comp  = trap.source(id_trap);
fprintf('Using trap depths: '); fprintf('%d ', z_comp); fprintf('m\n\n');

% ---------------------------------------------------------------
% 2. Run model (flux BC, best config)
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

k_bc   = 2;
dz     = col_grid.dz;
n_z    = col_grid.n_z;
dt     = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
k_plot = 2:10;
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

d_cm  = bc.d_model_um * 1e-4;
w_bin = (66 * d_cm .^ 0.62)';

sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);

% spinup
for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc, :) = Y(k_bc, :) + flux_src;
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        end
    end
    phi_after = mean(sum(Y + Yfp, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('Converged at cycle %d\n', icyc); break;
    end
end

% final run: accumulate BV flux on cast days
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);
F_bv_sum = zeros(n_z, 1);
n_cast   = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc, :) = Y(k_bc, :) + flux_src;
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        for k = 1:n_z
            F_bv_sum(k) = F_bv_sum(k) + sum(w_bin .* (Y(k,:) + Yfp(k,:)));
        end
        n_cast = n_cast + 1;
    end
end
F_bv = F_bv_sum / max(n_cast, 1);   % [m3 BV m-2 d-1] at each layer

% model BV flux at comparison depths
k_comp = zeros(1, numel(target_z));
for i = 1:numel(target_z)
    [~, k_comp(i)] = min(abs(col_grid.z_centers - z_comp(i)));
end
F_bv_comp = F_bv(k_comp);   % [m3 BV m-2 d-1]

% ---------------------------------------------------------------
% 3. Print table
% ---------------------------------------------------------------
fprintf('\n=== ESTAPA POC FLUX vs MODEL BV FLUX ===\n');
fprintf('%-8s  %-6s  %-16s  %-16s  %-12s  %-14s\n', ...
    'Depth', 'Src', 'Trap POC', 'Model F_bv', 'Ratio', 'Implied C:BV');
fprintf('%-8s  %-6s  %-16s  %-16s  %-12s  %-14s\n', ...
    '(m)', '', '(mg C m-2 d-1)', '(m3 m-2 d-1)', 'trap/model', '(mg C m-3 BV)');
for i = 1:numel(target_z)
    ratio = poc_comp(i) / max(F_bv_comp(i), 1e-30);
    fprintf('%5.0f m   %-6s  %12.1f      %12.2e      %8.2f    %12.2e\n', ...
        z_comp(i), src_comp{i}, poc_comp(i), F_bv_comp(i), ratio, ratio);
end

fprintf('\n-- POC flux profile (mg C m-2 d-1) --\n');
fprintf('All Estapa depths:\n');
for i = 1:numel(trap.depths)
    fprintf('  %5.0f m [%s]: %.1f +/- %.1f\n', ...
        trap.depths(i), trap.source{i}, trap.flux_POC(i), trap.flux_POC_sd(i));
end

% ---------------------------------------------------------------
% 4. Figure
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 14 14], 'Color', 'white');

subplot(1, 2, 1);
hold on;

% Estapa trap POC flux
errorbar(poc_comp, z_comp, zeros(size(poc_comp)), zeros(size(poc_comp)), ...
    poc_sd, poc_sd, 'bs', 'MarkerSize', 5, 'LineWidth', 1.2, ...
    'DisplayName', 'Estapa trap');

% model BV flux (right y-axis not available in simple plot — use separate panel)
% Here just show relative shape by normalizing model to match trap at 75m
F_mod_norm = F_bv(k_plot) * (poc_comp(1) / F_bv_comp(1));
z_mod = col_grid.z_centers(k_plot);
plot(F_mod_norm, z_mod, 'k-', 'LineWidth', 1.3, 'DisplayName', 'Model (scaled)');

set(gca, 'YDir', 'reverse', 'XScale', 'log', ...
    'YLim', [60 550], 'FontSize', 7);
xlabel('POC flux (mg C m^{-2} d^{-1})');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 6);
title('POC flux profile', 'FontWeight', 'normal');
hold off;

subplot(1, 2, 2);
% implied C:BV ratio at each comparison depth
ratio_vec = poc_comp ./ max(F_bv_comp, 1e-30);
hold on;
plot(ratio_vec, z_comp, 'bs-', 'MarkerSize', 5, 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse', 'XScale', 'log', ...
    'YLim', [60 550], 'FontSize', 7);
xlabel('Implied C:BV (mg C m^{-3} BV)');
ylabel('Depth (m)');
title('Implied POC per BV', 'FontWeight', 'normal');
hold off;

saveas(gcf, fullfile(fig_dir, 'estapa_comparison.png'));
fprintf('\nSaved estapa_comparison.png\n');

% ---------------------------------------------------------------
function cfg = cfg_best()
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

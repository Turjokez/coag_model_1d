% run_microbe_flux_sensitivity.m
% Test three microbe_r0 values: 0, 0.01, 0.02.
% Flux profiles normalized to 1 at 100 m (transfer efficiency).
%
% Shows that r0 ~ 0.014 is needed for realistic deep attenuation.
% Fitted value from inverse fit: r0=0.014 (100m BC), r0=0.020 (surface BC).

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));
addpath(fullfile(script_dir, '..', 'data'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

% --- paths ---
mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);
dz       = col_grid.dz;
z_c      = col_grid.z_centers;

% --- base config (best fit, alpha=0.093, Da*5, mining) ---
cfg_base = SimulationConfig();
cfg_base.n_sections     = 30;
cfg_base.sinking_law    = 'kriest_8';
cfg_base.ds_kernel_mode = 'sinking_law';
cfg_base.r_to_rg        = 1.6;
cfg_base.alpha          = 0.093;
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

% --- BC at 100 m, bc_scale=0.42 ---
bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, 2:10);
phi_bc_daily = bc.phi_bc_daily * 0.42;   % best fit bc_scale
n_days       = bc.n_days;

% --- r0 values to test ---
r0_vals = [0,    0.01,  0.02 ];
labels  = {'r_0 = 0 (no microbe)', 'r_0 = 0.01 d^{-1}', 'r_0 = 0.02 d^{-1}'};
colors  = {'b',  [0 0.55 0],  'r'};

[~, k100] = min(abs(z_c - 100));

F_store = zeros(col_grid.n_z, numel(r0_vals));

for ir = 1:numel(r0_vals)
    cfg = copy(cfg_base);
    if r0_vals(ir) > 0
        cfg.enable_microbe = true;
        cfg.microbe_r0     = r0_vals(ir);
    else
        cfg.enable_microbe = false;
    end

    sim   = ColumnSimulation(cfg, col_grid, prof);
    w_bin = 66 * sim.size_grid.dcomb(:)' .^ 0.62;

    dt            = 0.25;
    steps_per_day = round(1/dt);
    k_bc          = 2;
    spinup_tol    = 0.01;
    max_cycles    = 80;

    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);

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
            break
        end
    end
    fprintf('r0=%.3f: converged cycle %d\n', r0_vals(ir), icyc);

    F_bv = sum((Y + Yfp) .* w_bin, 2);
    F_store(:, ir) = F_bv / F_bv(k100);   % normalize to 1 at 100 m
end

% --- print TE table ---
fprintf('\n--- Transfer efficiency (relative to 100 m) ---\n');
fprintf('  %6s', 'z(m)');
for ir = 1:numel(r0_vals)
    fprintf('  %12s', sprintf('r0=%.2f', r0_vals(ir)));
end
fprintf('\n');
for zz = [200, 500, 975]
    [~, k] = min(abs(z_c - zz));
    fprintf('  %6.0f', z_c(k));
    for ir = 1:numel(r0_vals)
        fprintf('  %11.1f%%', F_store(k,ir)*100);
    end
    fprintf('\n');
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 10],'Color','white');
hold on;
for ir = 1:numel(r0_vals)
    plot(F_store(:,ir), z_c, '-', 'Color', colors{ir}, ...
        'LineWidth', 1.2, 'DisplayName', labels{ir});
end
xline(1,'k:','LineWidth',0.8,'HandleVisibility','off');
set(gca, 'YDir','reverse', 'XScale','log', 'FontSize',fs, ...
    'Box','on', 'YLim',[50 1000], 'XLim',[0.01 1.5]);
xlabel('flux / flux(100 m)', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','southeast', 'FontSize',fs, 'Box','off');

fig_dir  = fullfile(script_dir, '..', '..', 'docs', 'figs');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
exportgraphics(gcf, fullfile(fig_dir, 'microbe_flux_sensitivity.png'), 'Resolution',200);
fprintf('\nSaved microbe_flux_sensitivity.png\n');

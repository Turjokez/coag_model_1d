% run_shallow_flux_diagnosis.m
% Why does BV flux at 175 m exceed 100 m with r0=0?
%
% Isolate which process causes the shallow flux increase by
% turning off one physics term at a time.
%
% Runs:
%   A. full config (reference)
%   B. no zoo (enable_zoo=false)
%   C. no mining (enable_mining=false)
%   D. zoo but no fecal rerouting (zoo_p=0)
%   E. no zoo, no mining (minimal)
%
% All runs: r0=0, alpha=0.10, same BC.
% Metric: BV flux F(z) = sum(w .* Y) at each layer.
% If flux at 175m > flux at 100m in run E (minimal), it's sinking artifact.

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

% --- base config ---
cfg_base = SimulationConfig();
cfg_base.n_sections     = 30;
cfg_base.sinking_law    = 'kriest_8';
cfg_base.ds_kernel_mode = 'sinking_law';
cfg_base.r_to_rg        = 1.6;
cfg_base.alpha          = 0.10;
cfg_base.enable_coag    = true;
cfg_base.enable_disagg  = true;
cfg_base.disagg_mode    = 'operator_split';
cfg_base.disagg_dmax_A  = 9.39e-6 * 5;
cfg_base.enable_zoo     = true;
cfg_base.zoo_c          = 0.025;
cfg_base.zoo_s          = 1.3e-5;
cfg_base.zoo_p          = 0.5;
cfg_base.zoo_ic         = 7;
cfg_base.enable_microbe = false;
cfg_base.enable_mining  = true;

bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, 2:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

% Th-234 absolute scale at 100m
F_th234 = 4.65;
[~, k100] = min(abs(z_c - 100));

% --- define runs ---
run_names = {'A: full', 'B: no zoo', 'C: no mining', 'D: zoo, fp=0', 'E: no zoo/mining'};
n_runs = numel(run_names);
colors = {'b','r','g',[0.8 0.4 0],[0.5 0.5 0.5]};

F_store = zeros(col_grid.n_z, n_runs);

for ir = 1:n_runs
    cfg = copy(cfg_base);

    switch ir
        case 1  % A: full
            % no change
        case 2  % B: no zoo
            cfg.enable_zoo = false;
        case 3  % C: no mining
            cfg.enable_mining = false;
        case 4  % D: zoo but no fecal (fp removed, just grazing loss)
            cfg.zoo_p = 0;
        case 5  % E: minimal (no zoo, no mining)
            cfg.enable_zoo    = false;
            cfg.enable_mining = false;
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

    % BV flux at each layer, scaled to Th-234 at 100m
    F_bv = sum((Y + Yfp) .* w_bin, 2);
    scale = F_th234 / F_bv(k100);
    F_store(:, ir) = F_bv * scale;

    fprintf('%s: converged cycle %d, F(100m)=%.2f, F(175m)=%.2f, TE(175m)=%.0f%%\n', ...
        run_names{ir}, icyc, F_store(k100, ir), ...
        interp1(z_c, F_store(:,ir), 175, 'linear'), ...
        interp1(z_c, F_store(:,ir), 175, 'linear') / F_th234 * 100);
end

% --- figure: TE relative to 100m, shallow zone ---
fs = 7;
figure('Units','centimeters','Position',[2 2 10 12],'Color','white');
hold on;

% Martin reference
z_martin   = (100:5:400)';
TE_martin  = (z_martin / 100) .^ (-0.86);
plot(TE_martin, z_martin, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Martin b=0.86');

for ir = 1:n_runs
    TE = F_store(:, ir) / F_th234;
    plot(TE, z_c, '-', 'Color', colors{ir}, 'LineWidth', 1.2, ...
        'DisplayName', run_names{ir});
end

xline(1, 'k--', 'LineWidth', 0.5, 'HandleVisibility', 'off');
set(gca, 'YDir','reverse', 'XScale','log', 'FontSize',fs, ...
    'Box','on', 'YLim',[50 400], 'XLim',[0.3 3]);
xlabel('TE  (relative to 100 m)', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location','southeast', 'FontSize',fs, 'Box','off');
title('shallow flux diagnosis (100-400 m)', 'FontWeight','normal', 'FontSize',fs);

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figs');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
fig_path = fullfile(fig_dir, 'shallow_flux_diagnosis.png');
exportgraphics(gcf, fig_path, 'Resolution', 200);
fprintf('\nFigure saved: %s\n', fig_path);

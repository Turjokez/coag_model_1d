% run_example_03.m
% Example 3: full physics — zoo grazing, fecal pellets, mining.
%
% Builds on Example 2 by adding:
%   - Zooplankton grazing (Stemmann 2004 filter + flux feeders)
%   - Fecal pellet tracking (Yfp, separate sinking law)
%   - Micro-zooplankton mining
%
% Shows aggregate vs fecal pellet depth profiles on same axes.
%
% Saves: docs/figures/example_03_profile.png

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

% --- config: Example 2 + zoo + mining ---
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.r_to_rg        = 1.6;
cfg.alpha          = 0.10;
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_zoo     = true;       % <- filter + flux feeders
cfg.zoo_c          = 0.025;      % clearance rate
cfg.zoo_s          = 1.3e-5;     % capture cross-section
cfg.zoo_p          = 0.5;        % egestion fraction
cfg.zoo_ic         = 7;          % fecal bin -> bin 8 (~115 um)
cfg.enable_microbe = false;
cfg.enable_mining  = true;       % <- micro-zoo mining

sim   = ColumnSimulation(cfg, col_grid, prof);
d_cm  = sim.size_grid.dcomb(:)';
w_bin = 66 * d_cm .^ 0.62;

k_plot = 2:10;
bc = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;
spinup_tol    = 0.01;
max_cycles    = 80;

Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);

% spinup
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

% --- figure ---
fs   = 8;
xagg = sum(Y, 2);
xfp  = sum(Yfp, 2);
figure('Units','centimeters','Position',[2 2 14 8],'Color','white');

subplot(1,2,1);
plot(xagg, col_grid.z_centers, 'k-', 'LineWidth', 1.2);
set(gca, 'YDir','reverse', 'FontSize',fs, 'Box','on', 'YLim',[0 1000]);
xlabel('aggregate BV (m^3 m^{-3})', 'FontSize',fs);
ylabel('depth (m)', 'FontSize',fs);
title('(a) aggregates', 'FontWeight','normal','FontSize',fs);

subplot(1,2,2);
plot(xfp, col_grid.z_centers, 'k--', 'LineWidth', 1.2);
set(gca, 'YDir','reverse', 'FontSize',fs, 'Box','on', 'YLim',[0 1000], ...
    'YTickLabel',{});
xlabel('fecal BV (m^3 m^{-3})', 'FontSize',fs);
title('(b) fecal pellets', 'FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
print(gcf, fullfile(fig_dir,'example_03_profile.png'), '-dpng', '-r150');
fprintf('Saved example_03_profile.png\n');

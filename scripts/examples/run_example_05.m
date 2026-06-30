% run_example_05.m
% Example 5: full model comparison to EXPORTS-NA UVP.
%
% All physics on (coag + disagg + zoo + fecal + microbe + mining).
% Best-fit parameters: alpha=0.093, bc_scale=0.42, r0=0.014.
% Compare quasi-steady BV profile to UVP observations at 125, 325, 475 m.
%
% Saves: docs/figures/example_05_full_model.png

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));
addpath(fullfile(script_dir, '..', 'data'));
addpath(fullfile(script_dir, '..', 'inverse'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

% --- paths ---
mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% --- grid and profiles ---
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);
z_c      = col_grid.z_centers;
dz       = col_grid.dz;

% --- best-fit config (all physics on) ---
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.r_to_rg        = 1.6;
cfg.alpha          = 0.093;      % best fit June 2026
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_zoo     = true;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.enable_microbe = true;
cfg.microbe_r0     = 0.014;      % best fit June 2026
cfg.enable_mining  = true;

% --- BC at 100 m with bc_scale = 0.42 ---
bc_scale     = 0.42;             % best fit June 2026
bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 2:10);
phi_bc_daily = bc.phi_bc_daily * bc_scale;
n_days       = bc.n_days;

% --- spinup ---
Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);
spinup_tol = 0.01;
max_cycles  = 80;
dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;

sim   = ColumnSimulation(cfg, col_grid, prof);
d_cm  = sim.size_grid.dcomb(:)';
w_bin = 66 * d_cm .^ 0.62;
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
    phi_after  = mean(sum(Y + Yfp, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('Converged at cycle %d\n', icyc);
        break;
    end
end

% --- load UVP observations at comparison depths ---
obs_depths = [125, 325, 475];
obs        = load_uvp_obs(uvp_file, obs_depths);
bv_obs     = obs.bv_total;   % mean total BV [m3/m3] at each depth

% --- model BV at comparison depths ---
bv_mod = zeros(1, 3);
for k = 1:3
    [~, iz] = min(abs(z_c - obs_depths(k)));
    bv_mod(k) = sum(Y(iz,:) + Yfp(iz,:));
end

% --- print ratios ---
fprintf('\n');
for k = 1:3
    fprintf('model/obs at %dm: %.2f\n', obs_depths(k), bv_mod(k)/bv_obs(k));
end

% --- figure (2 panels) ---
bv_total = sum(Y + Yfp, 2);
fs = 8;

figure('Units','centimeters','Position',[2 2 14 9],'Color','white');

% left: BV profile vs depth (log x)
subplot(1,2,1);
pv = bv_total;  pv(pv <= 0) = NaN;
hold on;
plot(pv, z_c, 'k-', 'LineWidth', 1.4, 'DisplayName', 'model');
plot(bv_obs, obs_depths, 'ko', 'MarkerSize', 5, 'MarkerFaceColor','k', ...
    'DisplayName', 'UVP');
set(gca, 'YDir','reverse', 'XScale','log', 'FontSize',fs, 'Box','on', ...
    'YLim',[0 1000]);
xlabel('total BV (m^3 m^{-3})', 'FontSize',fs);
ylabel('depth (m)', 'FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('(a) BV profile', 'FontWeight','normal','FontSize',fs);

% right: model/obs ratio (horizontal bars)
subplot(1,2,2);
ratio = bv_mod ./ bv_obs;
barh(ratio, 0.5, 'FaceColor',[0.6 0.6 0.6], 'EdgeColor','k');
hold on;
plot([1 1], [0 4], 'k-', 'LineWidth', 0.8);
set(gca, 'YDir','reverse', 'YTick',1:3, ...
    'YTickLabel',{'125 m','325 m','475 m'}, 'FontSize',fs, 'Box','on', ...
    'XLim',[0 1.6], 'XTick',0:0.4:1.6);
xlabel('model / obs', 'FontSize',fs);
title('(b) model/obs', 'FontWeight','normal','FontSize',fs);

out_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
print(gcf, fullfile(out_dir,'example_05_full_model.png'), '-dpng', '-r150');
fprintf('Saved example_05_full_model.png\n');

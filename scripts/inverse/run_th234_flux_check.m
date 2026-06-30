% run_th234_flux_check.m
% Run 1-D column model with best config.
% Compute BV flux profile and scale to Th-234 benchmark at 100 m.
%
% Th-234 gives: F_POC(100 m) = 4.65 mmol/m^2/day (from bottle deficit).
% Model gives BV flux shape. We use Th-234 to set the absolute scale.
%
% ATI pump result (run_amaral_poc_inversion.m) gives 20.9 mmol/m^2/day at 95 m --
% 4.5x higher than Th-234 due to swimmer contamination in pump PL fraction.
%
% Config: alpha=0.10, Da*5, enable_mining=true, no microbe.

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

% --- best config (June 12 2026) ---
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
cfg.enable_zoo     = true;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.enable_microbe = false;
cfg.enable_mining  = true;

% --- boundary condition at 100 m from UVP ---
k_plot = 2:10;
bc           = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;

% --- set up model ---
sim   = ColumnSimulation(cfg, col_grid, prof);
w_bin = 66 * sim.size_grid.dcomb(:)' .^ 0.62;   % kriest_8 [m/day]

dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;      % insert flux at layer 2 (100 m BC)
spinup_tol    = 0.01;
max_cycles    = 80;

Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);

% --- spin up to quasi-steady state ---
fprintf('Spinning up...\n');
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
        fprintf('Converged at cycle %d\n', icyc);
        break
    end
end

% --- compute BV flux at each layer [m/day * m^3/m^3 = m^3/m^2/day] ---
F_bv = sum((Y + Yfp) .* w_bin, 2);
z_c  = col_grid.z_centers;

% --- scale to Th-234 at 100 m ---
% Th-234 bottle deficit gives F_POC(100 m) = 4.65 mmol/m^2/day.
% We use this as the absolute calibration for the model BV flux.
F_th234  = 4.65;     % mmol/m^2/day
z_th234  = 100;      % m

[~, k100] = min(abs(z_c - z_th234));
scale     = F_th234 / F_bv(k100);    % mmol / m^3 (implied POC per unit BV flux)
F_poc     = F_bv * scale;            % mmol/m^2/day (model, scaled to Th-234 at 100 m)

% --- ATI pump result (from report_june25_amaral_poc_inversion) ---
% These fluxes include swimmer contamination (4.5x Th-234 at 95 m).
wS_fit = 1.95; wL_fit = 16.8;
z_ati  = [95,    125,   175,   330,   500  ];
PS_ati = [0.945, 0.420, 0.285, 0.229, 0.233];
PL_ati = [1.133, 0.241, 0.115, 0.131, 0.069];
F_ati  = wS_fit * PS_ati + wL_fit * PL_ati;   % mmol/m^2/day

% --- Martin curve reference: F(z) = F_100 * (z/100)^(-0.86) ---
z_martin   = (100:10:1000)';
F_martin   = F_th234 * (z_martin / 100) .^ (-0.86);

% --- print summary ---
fprintf('\n--- Flux at key depths ---\n');
fprintf('  %6s  %10s  %10s  %6s\n', 'z(m)', 'model', 'Martin', 'ratio');
for zz = [100, 200, 500, 975]
    [~, k]   = min(abs(z_c   - zz));
    [~, km]  = min(abs(z_martin - zz));
    fprintf('  %6.0f  %10.2f  %10.2f  %6.1fx\n', ...
        z_c(k), F_poc(k), F_martin(km), F_poc(k)/F_martin(km));
end
fprintf('\n');
fprintf('  Th-234 at 100 m:  %.2f mmol/m^2/day\n', F_th234);
fprintf('  ATI pump at 95 m: %.2f mmol/m^2/day  (%.1fx Th-234 = swimmers)\n', ...
    F_ati(1), F_ati(1)/F_th234);
fprintf('  Model TE (975m / 100m): %.1f%%  (Martin: %.1f%%)\n', ...
    F_poc(end)/F_th234*100, F_martin(end)/F_th234*100);

% --- figure: 2 panels ---
fs = 7;
figure('Units','centimeters','Position',[2 2 16 10],'Color','white');

% left: absolute flux vs depth
subplot(1,2,1);
hold on;
plot(F_martin, z_martin, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Martin b=0.86');
plot(F_poc, z_c, 'b-', 'LineWidth', 1.2, 'DisplayName', 'model (Th-scaled)');
plot(F_ati, z_ati, 'r--s', 'MarkerSize', 4, 'LineWidth', 1.0, ...
    'DisplayName', 'ATI pump (+swimmers)');
plot(F_th234, z_th234, 'k^', 'MarkerSize', 7, 'MarkerFaceColor', 'k', ...
    'DisplayName', 'Th-234 (bottle)');
set(gca, 'YDir', 'reverse', 'XScale', 'log', 'FontSize', fs, ...
    'Box', 'on', 'YLim', [50 1000]);
xlabel('POC flux (mmol m^{-2} d^{-1})', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location', 'southeast', 'FontSize', fs, 'Box', 'off');

% right: TE relative to 100 m
subplot(1,2,2);
hold on;
TE_model  = F_poc / F_th234;
TE_ati    = F_ati / F_ati(1);
TE_martin = F_martin / F_th234;
plot(TE_martin, z_martin, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Martin b=0.86');
plot(TE_model, z_c, 'b-', 'LineWidth', 1.2, 'DisplayName', 'model');
plot(TE_ati,   z_ati, 'r--s', 'MarkerSize', 4, 'LineWidth', 1.0, ...
    'DisplayName', 'ATI pump');
set(gca, 'YDir', 'reverse', 'XScale', 'log', 'FontSize', fs, ...
    'Box', 'on', 'YLim', [50 1000]);
xlabel('TE  (relative to 100 m)', 'FontSize', fs);
ylabel('depth (m)', 'FontSize', fs);
legend('Location', 'southeast', 'FontSize', fs, 'Box', 'off');

% save
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figs');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
fig_path = fullfile(fig_dir, 'flux_th234_check.png');
exportgraphics(gcf, fig_path, 'Resolution', 200);
fprintf('\nFigure saved: %s\n', fig_path);

% run_data_column_daily.m
%
% Data-driven 1-D column run with spinup.
%
% Surface forcing: UVP phi from top 5 m, reset each day.
% Physics: real eps(z), T(z), S(z) from keps_for_dave.mat.
%
% Steps:
%   1. Build DepthProfile from keps data
%   2. Build daily surface phi from UVP
%   3. Spinup: repeat 26-day forcing until depth-mean phi converges (<1%)
%   4. Comparison run: one more 26-day pass from spun-up state
%   5. Plot: model depth profile vs UVP observed

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

% --- paths ---
mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% --- config (EXPORTS-ready: n=30, dt=0.25 day, 4 steps/day) ---
cfg = SimulationConfig();
cfg.n_sections      = 30;
cfg.sinking_law     = 'kriest_8';
cfg.disagg_mode     = 'operator_split';
cfg.disagg_dmax_cm  = 1.0;      % fallback; depth run uses eps(z)
cfg.enable_coag     = true;
cfg.enable_disagg   = true;
cfg.enable_zoo      = true;
cfg.enable_microbe  = true;
cfg.enable_mining   = true;
cfg.alpha           = 0.5;     % best-fit from 2D grid search
cfg.microbe_r0      = 0.03;    % best-fit from r0 scan (run_r0_scan.m)
cfg.surface_pp_mu   = 0.1;
cfg.r_to_rg         = 1.6;
cfg.zoo_c           = 0.025;    % Stemmann 2004
cfg.zoo_s           = 1.3e-5;   % Stemmann 2004
cfg.zoo_p           = 0.5;      % Stemmann 2004
cfg.zoo_ic          = 7;        % bin 8, about 115 um
cfg.mining_s        = 1.3e-5;
cfg.fp_alpha_cross  = 0.5;
cfg.validate();

dt          = 0.25;   % day (4 steps per day, CFL ~ 0.5)
steps_per_day = round(1 / dt);   % = 4

% --- column grid: 1000 m, 20 layers ---
col_grid = ColumnGrid(1000, 20);

% --- depth profile from real keps data ---
prof = load_keps(mat_path, col_grid.z_centers);
fprintf('DepthProfile built from keps data\n');
fprintf('  eps range: %.2e to %.2e cm^2/s^3\n', min(prof.eps), max(prof.eps));
fprintf('  T range: %.1f to %.1f C\n', min(prof.T_K)-273.15, max(prof.T_K)-273.15);

% --- daily surface phi from UVP ---
daily = get_daily_surface_phi(uvp_file, cfg, col_grid);
n_days = daily.n_days;
fprintf('Daily surface phi: %d days (%d with real UVP data)\n', ...
    n_days, sum(daily.has_data));

% --- set up simulation ---
sim = ColumnSimulation(cfg, col_grid, prof);

% check CFL with real sinking speeds
% w_z is m/day and dz is m, so CFL = w*dt/dz
w_max = max(sim.rhs.w_z(:));
cfl = ColumnTransport.maxCFL(sim.rhs.w_z, prof.Kz, col_grid.dz, dt);
fprintf('CFL check: w_max=%.2f m/day, dt=%.2f day, dz=%g m -> CFL=%.3f\n', ...
    w_max, dt, col_grid.dz, cfl);
if cfl > 0.9
    warning('CFL = %.3f > 0.9. Reduce dt.', cfl);
end

% --- spinup: repeat 26-day forcing until column reaches steady state ---
n_z   = col_grid.n_z;
n_sec = cfg.n_sections;

% start from zeros
Y   = zeros(n_z, n_sec);
Yfp = zeros(n_z, n_sec);

spinup_tol   = 0.01;   % converge when depth-mean phi changes < 1%
max_cycles   = 50;
conv_cycle   = nan;

fprintf('Spinup: repeating %d-day forcing until depth-mean phi < %.0f%% change...\n', ...
    n_days, spinup_tol*100);
t_start = tic;

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 3), 2);   % depth profile before this cycle

    for i_day = 1:n_days
        for i_step = 1:steps_per_day
            Y(1, :) = daily.phi(i_day, :);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(1, :) = daily.phi(i_day, :);
        end
    end

    phi_after = mean(sum(Y + Yfp, 3), 2);
    rel_change = max(abs(phi_after - phi_before) ./ max(phi_before, 1e-20));
    fprintf('  cycle %2d: max rel change = %.4f\n', icyc, rel_change);

    if rel_change < spinup_tol
        conv_cycle = icyc;
        fprintf('  converged at cycle %d.\n', icyc);
        break;
    end
end

if isnan(conv_cycle)
    fprintf('  warning: spinup did not converge in %d cycles.\n', max_cycles);
end

% --- comparison run from spun-up state ---
% storage: save every day
Y_daily   = zeros(n_days, n_z, n_sec);
Yfp_daily = zeros(n_days, n_z, n_sec);

fprintf('Comparison run: %d days from spun-up state...\n', n_days);

for i_day = 1:n_days
    for i_step = 1:steps_per_day
        Y(1, :) = daily.phi(i_day, :);
        [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
        Y(1, :) = daily.phi(i_day, :);
    end
    Y_daily(i_day, :, :)   = Y;
    Yfp_daily(i_day, :, :) = Yfp;
end

elapsed = toc(t_start);
fprintf('Done. Elapsed: %.1f s\n', elapsed);

% --- load UVP observed profile for comparison ---
uvp = parse_uvp(uvp_file);

% map UVP cruise-mean phi to model z grid
% UVP depth bins -> interpolate to model z
% use aggregate-sized UVP only; larger objects are mostly zooplankton
mask_agg = uvp.d_um >= 100 & uvp.d_um < 2000;
uvp_phi_clean = uvp.phi(:, mask_agg);
uvp_phi_clean(isnan(uvp_phi_clean)) = 0;
uvp_phi_total = sum(uvp_phi_clean, 2);   % n_uvp_depths x 1

% interpolate UVP phi to model z centers
uvp_phi_model = interp1(uvp.depth_m, uvp_phi_total, col_grid.z_centers, ...
    'pchip', 'extrap');
uvp_phi_model = max(0, uvp_phi_model);

% model cruise-mean phi (aggregate + fecal, sum over bins, mean over days)
model_phi_total = squeeze(sum(Y_daily + Yfp_daily, 3));   % n_days x n_z
model_phi_mean  = mean(model_phi_total, 1)';  % n_z x 1

% --- figures ---
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% Figure 1: model vs UVP depth profile (cruise mean)
% semilogx = log phi on x-axis, linear depth on y-axis
figure;
semilogx(model_phi_mean,    col_grid.z_centers, 'b-',  'DisplayName', 'model (daily forced)');
hold on;
semilogx(uvp_phi_model,     col_grid.z_centers, 'r--', 'DisplayName', 'UVP <2000 um');
set(gca, 'YDir', 'reverse');
xlabel('\phi_{total}  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('model vs UVP: cruise mean');
saveas(gcf, fullfile(fig_dir, 'data_daily_depth_profile.png'));

% Figure 2: surface phi time series (model top layer vs UVP forcing)
model_surf = squeeze(Y_daily(:, 1, :) + Yfp_daily(:, 1, :));   % n_days x n_sec
model_surf_total = sum(model_surf, 2);

figure;
plot(daily.day_num, sum(daily.phi, 2), 'r-o', 'MarkerSize', 3, ...
    'DisplayName', 'UVP surface (forcing)');
hold on;
plot(daily.day_num, model_surf_total, 'b-', 'DisplayName', 'model surface (end of day)');
xlabel('day');
ylabel('\phi  [cm^3 cm^{-3}]');
legend;
title('surface phi: forcing vs model');
saveas(gcf, fullfile(fig_dir, 'data_daily_surface_time.png'));

% Figure 3: time-depth of model total phi (log10 scale)
% Shows how the column evolves over the cruise.
model_phi_log = log10(max(model_phi_total, 1e-12));   % n_days x n_z
figure;
imagesc(daily.day_num, col_grid.z_centers, model_phi_log');
set(gca, 'YDir', 'reverse');
xlabel('cruise day');
ylabel('depth  [m]');
colorbar;
title('model \phi (log_{10} cm^3 cm^{-3})');
saveas(gcf, fullfile(fig_dir, 'data_daily_timedepth.png'));

% Figure 4: selected-day profiles vs UVP cruise mean
% Pick early, middle, and late in the cruise.
sel_idx = round([1, n_days/2, n_days]);
sel_styles = {'b-', 'g-', 'r-'};
figure;
for k = 1:3
    id = sel_idx(k);
    semilogx(model_phi_total(id,:)', col_grid.z_centers, sel_styles{k}, ...
        'DisplayName', sprintf('model day %d', daily.day_num(id)));
    hold on;
end
semilogx(uvp_phi_model, col_grid.z_centers, 'k--', 'DisplayName', 'UVP mean <2000 \mum');
set(gca, 'YDir', 'reverse');
xlabel('\phi  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('model: day 1 / mid / end vs UVP mean');
saveas(gcf, fullfile(fig_dir, 'data_daily_profiles_selected.png'));

% --- summary ---
fprintf('\n--- Summary ---\n');
fprintf('UVP surface phi (mean):   %.3e cm^3/cm^3\n', mean(sum(daily.phi, 2)));
fprintf('Model surface phi (mean): %.3e cm^3/cm^3\n', mean(model_surf_total));
check_depths = [25, 75, 175, 300, 500, 975];
for zd = check_depths
    [~, iz] = min(abs(col_grid.z_centers - zd));
    ratio = model_phi_mean(iz) / max(uvp_phi_model(iz), 1e-20);
    fprintf('  z=%4.0f m  model=%.2e  UVP=%.2e  ratio=%.2f\n', ...
        col_grid.z_centers(iz), model_phi_mean(iz), uvp_phi_model(iz), ratio);
end

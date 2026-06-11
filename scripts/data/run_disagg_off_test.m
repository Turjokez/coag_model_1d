% run_disagg_off_test.m
%
% Test: does turning off disaggregation create more large particles?
%
% Runs two versions of the model:
%   (A) full physics (same as run_compare_spectrum)
%   (B) disaggregation OFF (enable_disagg = false)
%
% Compares the size spectrum at selected depths for the same cast day.
% Goal: if disagg-off creates large particles -> disagg is destroying them.
%       if disagg-off still has no large particles -> coagulation is the problem.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- base config (same as run_compare_spectrum) ---
cfg_base = SimulationConfig();
cfg_base.n_sections    = 30;
cfg_base.sinking_law   = 'kriest_8';
cfg_base.disagg_mode   = 'operator_split';
cfg_base.disagg_dmax_cm = 1.0;
cfg_base.enable_coag   = true;
cfg_base.enable_disagg = true;    % will override below for case B
cfg_base.enable_zoo    = true;
cfg_base.enable_microbe = true;
cfg_base.enable_mining  = true;
cfg_base.alpha          = 0.5;
cfg_base.microbe_r0     = 0.03;
cfg_base.surface_pp_mu  = 0.1;
cfg_base.r_to_rg        = 1.6;
cfg_base.zoo_c          = 0.025;
cfg_base.zoo_s          = 1.3e-5;
cfg_base.zoo_p          = 0.5;
cfg_base.zoo_ic         = 7;
cfg_base.mining_s       = 1.3e-5;
cfg_base.fp_alpha_cross = 0.5;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

col_grid = ColumnGrid(1000, 20);
prof     = load_keps(mat_path, col_grid.z_centers);
daily    = get_daily_surface_phi(uvp_file, cfg_base, col_grid);
n_days   = daily.n_days;
n_z      = col_grid.n_z;

% model bin diameters and widths [mm]
grid_cfg   = cfg_base.derive();
r_cm       = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_mm = (2 * r_cm * 1e4)' / 1000;   % 1 x n_sec, bin centers [mm]

% bin boundaries (geometric mean of adjacent centers)
log_d      = log(d_model_mm);
log_bnd    = [log_d(1)-(log_d(2)-log_d(1))/2, ...
              (log_d(1:end-1)+log_d(2:end))/2, ...
              log_d(end)+(log_d(end)-log_d(end-1))/2];
d_model_bounds_mm = exp(log_bnd);              % 1 x (n_sec+1)
dw_model_mm = diff(d_model_bounds_mm);         % 1 x n_sec, bin widths [mm]

% run one case: full config or disagg-off
function [model_phi] = run_case(cfg, col_grid, prof, daily, dt, steps_per_day, spinup_tol, max_cycles)
    n_z   = col_grid.n_z;
    n_sec = cfg.n_sections;
    n_days = daily.n_days;
    cfg.validate();
    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(n_z, n_sec);
    Yfp = zeros(n_z, n_sec);
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 3), 2);
        for i_day = 1:n_days
            for i_step = 1:steps_per_day
                Y(1, :) = daily.phi(i_day, :);
                [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
                Y(1, :) = daily.phi(i_day, :);
            end
        end
        phi_after  = mean(sum(Y + Yfp, 3), 2);
        rel_change = max(abs(phi_after - phi_before) ./ max(phi_before, 1e-20));
        if rel_change < spinup_tol, break; end
    end
    Y_daily   = zeros(n_days, n_z, n_sec);
    Yfp_daily = zeros(n_days, n_z, n_sec);
    for i_day = 1:n_days
        for i_step = 1:steps_per_day
            Y(1, :) = daily.phi(i_day, :);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(1, :) = daily.phi(i_day, :);
        end
        Y_daily(i_day, :, :)   = Y;
        Yfp_daily(i_day, :, :) = Yfp;
    end
    model_phi = Y_daily + Yfp_daily;   % n_days x n_z x n_sec
end

% --- case A: full physics ---
fprintf('Case A: full physics...\n');
cfg_A = cfg_base;
phi_A = run_case(cfg_A, col_grid, prof, daily, dt, steps_per_day, spinup_tol, max_cycles);
fprintf('  done.\n');

% --- case B: disagg off ---
fprintf('Case B: disagg off...\n');
cfg_B = cfg_base;
cfg_B.enable_disagg = false;
phi_B = run_case(cfg_B, col_grid, prof, daily, dt, steps_per_day, spinup_tol, max_cycles);
fprintf('  done.\n');

% --- parse UVP by date/depth for comparison ---
uvpd = parse_uvp_daily(uvp_file);
uvp_bin_mask = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_mm     = uvpd.d_um(uvp_bin_mask) / 1000;
dw_filt      = uvpd.dw(uvp_bin_mask);

% pick the best cast day (highest surface forcing)
[cast_dates_matched, ia, ib] = intersect(daily.dates, uvpd.dates);
[~, best] = max(sum(daily.phi(ia, :), 2));
id_A   = ia(best);
id_uvp = ib(best);
date_str = num2str(cast_dates_matched(best));
lbl = [date_str(5:6) '-' date_str(7:8)];

% check depths
check_z = [75, 200, 400];
styles_A   = {'b-',  'r-',  'g-'};
styles_B   = {'b--', 'r--', 'g--'};
styles_uvp = {'b:', 'r:', 'g:'};

% --- Figure 1: spectrum at three depths, full vs disagg-off vs UVP ---
figure;
hold on;
for kd = 1:3
    [~, iz_mod] = min(abs(col_grid.z_centers - check_z(kd)));
    [~, iz_uvp] = min(abs(uvpd.depth_m - check_z(kd)));

    % UVP spectrum [ppmV/mm] — ensure row vectors for consistent indexing
    phi_u = squeeze(uvpd.phi(id_uvp, iz_uvp, uvp_bin_mask));
    S_uvp = phi_u(:)' ./ dw_filt(:)' * 1e9;   % 1 x n_uvp_bins
    S_uvp(S_uvp <= 0) = NaN;

    % case A: full physics -> spectrum [ppmV/mm] = phi / dw [mm] * 1e9
    phi_a = squeeze(phi_A(id_A, iz_mod, :))';   % 1 x n_sec
    S_a   = phi_a ./ dw_model_mm * 1e6;

    % case B: disagg off
    phi_b = squeeze(phi_B(id_A, iz_mod, :))';
    S_b   = phi_b ./ dw_model_mm * 1e6;

    ok_a = S_a > 0 & isfinite(S_a);
    ok_b = S_b > 0 & isfinite(S_b);
    ok_u = ~isnan(S_uvp) & S_uvp > 0;

    if any(ok_a)
        loglog(d_model_mm(ok_a), S_a(ok_a), styles_A{kd}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('full  %dm', check_z(kd)));
    end
    if any(ok_b)
        loglog(d_model_mm(ok_b), S_b(ok_b), styles_B{kd}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('no disagg  %dm', check_z(kd)));
    end
    if any(ok_u)
        loglog(d_uvp_mm(ok_u), S_uvp(ok_u), styles_uvp{kd}, 'LineWidth', 1.5, ...
            'DisplayName', sprintf('UVP  %dm', check_z(kd)));
    end
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([0.05 10]);
ylim([1e-2 1e3]);
xlabel('ESD (mm)');
ylabel('ppmV mm^{-1}');
legend('location', 'southwest', 'FontSize', 6);
title(sprintf('disagg-off test  date=%s  (solid=full, dash=no disagg, dot=UVP)', lbl));
saveas(gcf, fullfile(fig_dir, 'disagg_off_spectrum.png'));

% --- Figure 2: cruise-mean depth profiles, full vs disagg-off ---
phi_A_mean = squeeze(mean(sum(phi_A, 3), 1))';   % n_z x 1
phi_B_mean = squeeze(mean(sum(phi_B, 3), 1))';

figure;
semilogx(phi_A_mean, col_grid.z_centers, 'b-',  'DisplayName', 'full physics');
hold on;
semilogx(phi_B_mean, col_grid.z_centers, 'r--', 'DisplayName', 'disagg off');
set(gca, 'YDir', 'reverse');
xlabel('\phi_{total}  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('disagg-off: total phi vs depth');
saveas(gcf, fullfile(fig_dir, 'disagg_off_profile.png'));

fprintf('Figures saved to %s\n', fig_dir);
fprintf('\nResult interpretation:\n');
fprintf('  If no-disagg has MORE large particles -> disagg is destroying large aggregates.\n');
fprintf('  If no-disagg still lacks large particles -> coagulation is not building them.\n');

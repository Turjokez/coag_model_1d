% run_dmax_sensitivity.m
%
% Sweep Dmax_A x [1, 3, 5, 10] and compare size spectrum vs UVP.
% D_max = Dmax_A * eps^(-1/4).  Default Dmax_A = 9.39e-6 m.
% At surface eps ~1e-6 m^2/s^3: Dmax_A x1 -> D_max~0.30 mm,
%   x3 -> 0.89 mm, x5 -> 1.49 mm, x10 -> 2.97 mm.
%
% Plot: loglog spectrum at z=75m and z=200m vs UVP best-cast day.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Base config
% ---------------------------------------------------------------
cfg = SimulationConfig();
cfg.n_sections      = 30;
cfg.sinking_law     = 'kriest_8';
cfg.disagg_mode     = 'operator_split';
cfg.disagg_dmax_cm  = 1.0;
cfg.enable_coag     = true;
cfg.enable_disagg   = true;
cfg.enable_zoo      = true;
cfg.enable_microbe  = true;
cfg.enable_mining   = true;
cfg.alpha           = 0.5;
cfg.microbe_r0      = 0.03;
cfg.microbe_use_temp = true;
cfg.microbe_tref_C  = 20;
cfg.surface_pp_mu   = 0.1;
cfg.r_to_rg         = 1.6;
cfg.zoo_c           = 0.025;
cfg.zoo_s           = 1.3e-5;
cfg.zoo_p           = 0.5;
cfg.zoo_ic          = 7;
cfg.mining_s        = 1.3e-5;
cfg.fp_alpha_cross  = 0.5;
cfg.validate();

col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);   % struct: .eps, .dates
daily     = get_daily_surface_phi(uvp_file, cfg, col_grid);
n_days    = daily.n_days;
n_z       = col_grid.n_z;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;
check_z       = [75, 200];   % depths for comparison

% ---------------------------------------------------------------
% 2. Model bin geometry (mm)
% ---------------------------------------------------------------
grid_cfg    = cfg.derive();
r_cm        = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_mm  = (2 * r_cm * 1e4)' / 1000;   % 1 x n_sec, mm
log_d       = log(d_model_mm);
log_bnd     = [log_d(1)-(log_d(2)-log_d(1))/2, ...
               (log_d(1:end-1)+log_d(2:end))/2, ...
               log_d(end)+(log_d(end)-log_d(end-1))/2];
dw_model_mm = diff(exp(log_bnd));           % 1 x n_sec, mm

% ---------------------------------------------------------------
% 3. UVP parse and filter (100-2000 um)
% ---------------------------------------------------------------
uvpd         = parse_uvp_daily(uvp_file);
uvp_bin_mask = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_mm     = uvpd.d_um(uvp_bin_mask) / 1000;
dw_filt_um   = uvpd.dw(uvp_bin_mask);      % um

% best cast day (highest surface phi that has a UVP cast)
[~, ia, ib]  = intersect(daily.dates, uvpd.dates);
[~, best]    = max(sum(daily.phi(ia, :), 2));
id_model     = ia(best);
id_uvp       = ib(best);

% UVP spectrum at each check depth
S_uvp = NaN(numel(check_z), sum(uvp_bin_mask));
for ip = 1:numel(check_z)
    [~, iz_u] = min(abs(uvpd.depth_m - check_z(ip)));
    phi_u = squeeze(uvpd.phi(id_uvp, iz_u, uvp_bin_mask));
    S_uvp(ip, :) = phi_u(:)' ./ dw_filt_um(:)' * 1e9;   % ppmV/mm (dw in um)
end

% ---------------------------------------------------------------
% 4. Dmax_A base value and scale setup
% ---------------------------------------------------------------
prof        = load_keps(mat_path, col_grid.z_centers);
Dmax_A_base = 9.39e-6;   % m, Parker default

scales  = [1, 3, 5, 10];
colors  = {'k', 'b', 'r', 'm'};
n_sc    = numel(scales);

% D_max at surface eps (1e-6 m^2/s^3) for legend labels
eps_surf_m2s3 = 1e-6;   % m^2/s^3, typical mixed layer
Dmax_mm  = Dmax_A_base * scales * eps_surf_m2s3^(-0.25) * 1000;  % mm

% store spectra: n_sc x n_depths x n_sec
S_mod = NaN(n_sc, numel(check_z), cfg.n_sections);

% ---------------------------------------------------------------
% 5. Run model once per scale
% ---------------------------------------------------------------
for ks = 1:n_sc
    sc = scales(ks);
    fprintf('Dmax_A x%d  (D_max at surface = %.2f mm) ...\n', sc, Dmax_mm(ks));

    cfg.disagg_dmax_A = Dmax_A_base * sc;
    cfg.validate();
    sim = ColumnSimulation(cfg, col_grid, prof);

    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            for i_step = 1:steps_per_day
                Y(1,:) = daily.phi(i_day,:);
                [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
                Y(1,:) = daily.phi(i_day,:);
            end
        end
        phi_after  = mean(sum(Y + Yfp, 2));
        rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
        if rel_change < spinup_tol
            fprintf('  converged at cycle %d\n', icyc);
            break;
        end
    end

    % one more pass to capture best-cast day snapshot
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    phi_snap = zeros(n_z, cfg.n_sections);
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(1,:) = daily.phi(i_day,:);
            [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
            Y(1,:) = daily.phi(i_day,:);
        end
        if i_day == id_model
            phi_snap = Y + Yfp;   % n_z x n_sec
        end
    end

    % extract spectrum at each check depth
    for ip = 1:numel(check_z)
        [~, iz_m] = min(abs(col_grid.z_centers - check_z(ip)));
        phi_row = reshape(phi_snap(iz_m, :), 1, []);
        S_mod(ks, ip, :) = phi_row ./ dw_model_mm * 1e6;   % ppmV/mm (dw in mm)
    end
end

% ---------------------------------------------------------------
% 6. Plot: 1 x 2 subplots (75m, 200m)
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 18 8]);

for ip = 1:numel(check_z)
    subplot(1, 2, ip);
    hold on;

    % UVP
    S_u = S_uvp(ip, :);
    S_u(S_u <= 0) = NaN;
    ok  = ~isnan(S_u);
    loglog(d_uvp_mm(ok), S_u(ok), 'k--', 'LineWidth', 1.5, 'DisplayName', 'UVP');

    % model cases
    for ks = 1:n_sc
        S_m = reshape(S_mod(ks, ip, :), 1, []);
        ok_m = S_m > 0 & isfinite(S_m);
        if any(ok_m)
            loglog(d_model_mm(ok_m), S_m(ok_m), colors{ks}, 'LineWidth', 1.2, ...
                'DisplayName', sprintf('x%d (%.2f mm)', scales(ks), Dmax_mm(ks)));
        end
    end

    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([0.05 10]);
    ylim([1e-2 1e3]);
    xlabel('ESD (mm)');
    ylabel('ppmV mm^{-1}');
    legend('location', 'southwest', 'FontSize', 6);
    title(sprintf('z = %d m', check_z(ip)));
end

saveas(gcf, fullfile(fig_dir, 'dmax_sensitivity.png'));
fprintf('Saved dmax_sensitivity.png\n');

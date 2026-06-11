% run_r0_scan.m
%
% Quick 1-D scan over microbe_r0 to find value that minimizes
% model/UVP ratio across the water column.
%
% Fixes: alpha = 0.5 (best-fit from 2D grid search)
% Scans: r0 in r0_vals
% Loss:  mean of squared log10(model/UVP) over check depths

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

r0_vals = [0.005, 0.01, 0.015, 0.02, 0.03, 0.05];

% column setup (same as run_data_column_daily)
col_grid = ColumnGrid(1000, 20);
prof     = load_keps(mat_path, col_grid.z_centers);
dt       = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

% UVP target (same filter as main script)
uvp = parse_uvp(uvp_file);
mask_agg = uvp.d_um >= 100 & uvp.d_um < 2000;
uvp_phi_clean = uvp.phi(:, mask_agg);
uvp_phi_clean(isnan(uvp_phi_clean)) = 0;
uvp_phi_total = sum(uvp_phi_clean, 2);
% depths to check (indices into col_grid.z_centers)
check_z = [75, 175, 275, 475, 750];
uvp_at_check = interp1(uvp.depth_m, uvp_phi_total, check_z, 'pchip', 'extrap');
uvp_at_check = max(0, uvp_at_check);

fprintf('%-8s', 'r0');
for k = 1:length(check_z)
    fprintf('  z=%3dm', check_z(k));
end
fprintf('    loss\n');

losses = nan(size(r0_vals));
for ir = 1:length(r0_vals)
    r0 = r0_vals(ir);

    % build config
    cfg = SimulationConfig();
    cfg.n_sections   = 30;
    cfg.sinking_law  = 'kriest_8';
    cfg.disagg_mode  = 'operator_split';
    cfg.disagg_dmax_cm = 1.0;
    cfg.enable_coag  = true;
    cfg.enable_disagg = true;
    cfg.enable_zoo   = true;
    cfg.enable_microbe = true;
    cfg.enable_mining  = true;
    cfg.alpha        = 0.5;
    cfg.microbe_r0   = r0;
    cfg.surface_pp_mu = 0.1;
    cfg.r_to_rg      = 1.6;
    cfg.zoo_c        = 0.025;
    cfg.zoo_s        = 1.3e-5;
    cfg.zoo_p        = 0.5;
    cfg.zoo_ic       = 7;
    cfg.mining_s     = 1.3e-5;
    cfg.fp_alpha_cross = 0.5;
    cfg.validate();

    daily = get_daily_surface_phi(uvp_file, cfg, col_grid);
    n_days = daily.n_days;
    n_z   = col_grid.n_z;
    n_sec = cfg.n_sections;

    sim = ColumnSimulation(cfg, col_grid, prof);

    % spinup
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
        phi_after = mean(sum(Y + Yfp, 3), 2);
        rel_change = max(abs(phi_after - phi_before) ./ max(phi_before, 1e-20));
        if rel_change < spinup_tol
            break;
        end
    end

    % one comparison pass
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

    model_phi_mean = mean(squeeze(sum(Y_daily + Yfp_daily, 3)), 1)';  % n_z x 1
    model_at_check = interp1(col_grid.z_centers, model_phi_mean, check_z, 'pchip', 'extrap');

    % loss = mean squared log10 ratio at check depths
    ratios = model_at_check ./ max(uvp_at_check, 1e-20);
    loss   = mean((log10(ratios)).^2);
    losses(ir) = loss;

    fprintf('%-8.4f', r0);
    for k = 1:length(check_z)
        fprintf('  %5.2f ', ratios(k));
    end
    fprintf('    %.4f\n', loss);
end

[~, ibest] = min(losses);
fprintf('\nBest r0 = %.4f  (loss = %.4f)\n', r0_vals(ibest), losses(ibest));

% run_mass_budget_check.m
%
% Find where model mass differs from UVP.
%
% The key question (Adrian, June 11): model gives 30-40x more mass than
% UVP at depth when running without the overlap fix. Find out why.
%
% Strategy: after spinup, compare total phi at each depth for:
%   (a) model -- all 30 bins (full size range 1 um to ~10 mm)
%   (b) model -- UVP range only (100-2000 um)
%   (c) UVP measured phi (100-2000 um, from cast data)
%
% If (b)/(c) ~ 1  -> the mismatch was only because we compared different
%                    size ranges (model all-bins vs UVP-range).
% If (b)/(c) >> 1 -> real mass overproduction in the UVP size range.
%
% Also prints: fraction of model mass in small bins (<100 um) per depth.
% That tells us if the small bins are the hidden source.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% ---------------------------------------------------------------
% 1. Config (same as run_compare_spectrum)
% ---------------------------------------------------------------
cfg = SimulationConfig();
cfg.n_sections       = 30;
cfg.sinking_law      = 'kriest_8';
cfg.disagg_mode      = 'operator_split';
cfg.disagg_dmax_cm   = 1.0;
cfg.disagg_dmax_A    = 9.39e-6 * 5;   % x5 as used in comparison run
cfg.enable_coag      = true;
cfg.enable_disagg    = true;
cfg.enable_zoo       = true;
cfg.enable_microbe   = true;
cfg.enable_mining    = true;
cfg.alpha            = 0.5;
cfg.microbe_r0       = 0.03;
cfg.microbe_use_temp = true;
cfg.microbe_tref_C   = 20;
cfg.surface_pp_mu    = 0.1;
cfg.r_to_rg          = 1.6;
cfg.zoo_c            = 0.025;
cfg.zoo_s            = 1.3e-5;
cfg.zoo_p            = 0.5;
cfg.zoo_ic           = 7;
cfg.mining_s         = 1.3e-5;
cfg.fp_alpha_cross   = 0.5;
cfg.validate();

col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);
daily     = get_daily_surface_phi(uvp_file, cfg, col_grid);
n_days    = daily.n_days;
n_z       = col_grid.n_z;

% ---------------------------------------------------------------
% 2. Model bin geometry
% ---------------------------------------------------------------
grid_cfg   = cfg.derive();
r_cm       = (0.75 / pi * grid_cfg.av_vol(:)).^(1/3);
d_model_um = 2 * r_cm * 1e4;   % diameter in um, n_sec x 1

% bin masks
mask_uvp  = d_model_um >= 100 & d_model_um <= 2000;   % UVP range
mask_small = d_model_um < 100;                          % below UVP detection
fprintf('Model bins in UVP range (100-2000 um): %d of %d\n', ...
    sum(mask_uvp), cfg.n_sections);
fprintf('Model bins below 100 um:               %d of %d\n', ...
    sum(mask_small), cfg.n_sections);

% ---------------------------------------------------------------
% 3. UVP cast data at each depth (from parse_uvp_daily)
% ---------------------------------------------------------------
uvpd         = parse_uvp_daily(uvp_file);
uvp_mask     = uvpd.d_um >= 100 & uvpd.d_um < 2000;
dw_uvp_um    = uvpd.dw(uvp_mask);   % bin widths, um

% find best cast day (highest surface phi with a UVP cast)
[~, ia, ib] = intersect(daily.dates, uvpd.dates);
[~, best]   = max(sum(daily.phi(ia, :), 2));
id_model    = ia(best);
id_uvp      = ib(best);
fprintf('\nBest cast day: %d (index %d in model, %d in UVP)\n', ...
    daily.dates(id_model), id_model, id_uvp);

% UVP phi at each depth for this cast (cm3/cm3, summed across UVP range)
uvp_phi_depth = zeros(n_z, 1);
for k = 1:n_z
    [~, iz_u] = min(abs(uvpd.depth_m - col_grid.z_centers(k)));
    phi_uvp_bins = squeeze(uvpd.phi(id_uvp, iz_u, uvp_mask));
    uvp_phi_depth(k) = sum(phi_uvp_bins(:));
end

% ---------------------------------------------------------------
% 4. Run model to steady state (spinup)
% ---------------------------------------------------------------
dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);

fprintf('Running spinup...\n');
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
        fprintf('Converged at cycle %d\n', icyc);
        break;
    end
end

% one more pass: capture snapshot on best cast day
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);
Y_snap = [];
for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    for i_step = 1:steps_per_day
        Y(1,:) = daily.phi(i_day,:);
        [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
        Y(1,:) = daily.phi(i_day,:);
    end
    if i_day == id_model
        Y_snap   = Y + Yfp;   % total (agg + fecal)
        Y_agg    = Y;
        Y_fp     = Yfp;
    end
end

% ---------------------------------------------------------------
% 5. Budget table per depth
% ---------------------------------------------------------------
fprintf('\n--- Mass Budget by Depth (best cast day: %d) ---\n', daily.dates(id_model));
fprintf('%-8s  %-12s  %-12s  %-12s  %-8s  %-8s\n', ...
    'Depth(m)', 'Model-All', 'Model-UVP', 'UVP', 'Ratio', 'SmallFrac');
fprintf('%-8s  %-12s  %-12s  %-12s  %-8s  %-8s\n', ...
    '', '[cm3/cm3]', '[cm3/cm3]', '[cm3/cm3]', 'Mod/UVP', '<100um');
fprintf('%s\n', repmat('-', 1, 68));

for k = 1:n_z
    mod_all   = sum(Y_snap(k,:));           % all 30 bins
    mod_uvp   = sum(Y_snap(k, mask_uvp));   % only UVP-range bins
    mod_small = sum(Y_snap(k, mask_small)); % below 100 um bins
    uvp_val   = uvp_phi_depth(k);

    if uvp_val > 0
        ratio = mod_uvp / uvp_val;
    else
        ratio = NaN;
    end
    if mod_all > 0
        small_frac = mod_small / mod_all;
    else
        small_frac = 0;
    end

    fprintf('%8.1f  %12.3e  %12.3e  %12.3e  %8.2f  %8.1f%%\n', ...
        col_grid.z_centers(k), mod_all, mod_uvp, uvp_val, ratio, 100*small_frac);
end

% ---------------------------------------------------------------
% 6. Summary statistics
% ---------------------------------------------------------------
fprintf('\n--- Summary ---\n');
mod_uvp_col  = sum(Y_snap(:, mask_uvp),  2);
uvp_col      = uvp_phi_depth;
ok           = uvp_col > 0 & mod_uvp_col > 0;
ratio_vec    = mod_uvp_col(ok) ./ uvp_col(ok);
fprintf('Median ratio (model UVP-range / UVP):  %.1fx\n', median(ratio_vec));
fprintf('Mean   ratio (model UVP-range / UVP):  %.1fx\n', mean(ratio_vec));

small_frac_mean = mean(sum(Y_snap(:, mask_small), 2) ./ max(sum(Y_snap, 2), 1e-30));
fprintf('Mean fraction of model mass in <100 um bins: %.1f%%\n', 100*small_frac_mean);

% surface BC check: how much mass do we set at surface per day?
surf_phi_total    = sum(daily.phi(id_model,:));
surf_phi_uvp_range = sum(daily.phi(id_model, mask_uvp));
surf_phi_small    = sum(daily.phi(id_model, mask_small));
fprintf('\nSurface BC on best cast day:\n');
fprintf('  Total phi (all bins):   %.3e cm3/cm3\n', surf_phi_total);
fprintf('  UVP-range bins only:    %.3e cm3/cm3  (%.1f%% of total)\n', ...
    surf_phi_uvp_range, 100*surf_phi_uvp_range/max(surf_phi_total,1e-30));
fprintf('  Small bins (<100 um):   %.3e cm3/cm3  (%.1f%% of total)\n', ...
    surf_phi_small, 100*surf_phi_small/max(surf_phi_total,1e-30));
fprintf('\nNote: small bins at surface come from get_daily_surface_phi bin-mapping.\n');
fprintf('If they are zero, model has no small particles at surface -> \n');
fprintf('all sub-100um mass at depth is from disaggregation of large particles.\n');

% ---------------------------------------------------------------
% 7. Simple figure: model vs UVP total phi profile
% ---------------------------------------------------------------
figure('Units','centimeters','Position',[2 2 10 14]);
hold on;
plot(mod_uvp_col * 1e6, col_grid.z_centers, 'b-', 'LineWidth', 1.5, ...
    'DisplayName', 'model (UVP range)');
plot(uvp_col * 1e6, col_grid.z_centers, 'k--', 'LineWidth', 1.5, ...
    'DisplayName', 'UVP measured');
set(gca, 'YDir', 'reverse');
xlabel('\phi [ppmV]');
ylabel('Depth (m)');
legend('location', 'southeast');
title('Total phi: model vs UVP');

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'mass_budget_profile.png'));
fprintf('\nSaved mass_budget_profile.png\n');

% run_bianchi_lite_test.m
%
% Test Bianchi-lite DVM: gut pool + body pool mortality source into Y.
%
% Three cases:
%   1. baseline (no DVM)
%   2. gut-only DVM
%   3. Bianchi-lite DVM (gut + body mortality to Y)
%
% Main check:
%   Does the deep 375-475 m ratio improve when B_mig adds POC to Y?

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
% 1. Setup
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

k_bc   = 2;
k_plot = 2:10;
z_mod  = col_grid.z_centers(k_plot);

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

bc           = get_daily_bc_at_depth(uvp_file, cfg_base('off'), col_grid, 100, 3:10);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% UVP reference: mean BV over cast days, 100-2000 um
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ia, ib] = intersect(bc.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    for ki = 1:numel(k_plot)
        z_target = z_mod(ki);
        [~, iz] = min(abs(uvpd.depth_m(mask_z) - z_target));
        phi_u = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
        if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz, :));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);

% Spectrum at 375 m
cfg_tmp  = cfg_base('off');
sg       = DerivedGrid(cfg_tmp);
d_um_mod = sg.dcomb * 1e4;
mask_mod = d_um_mod >= 100 & d_um_mod < 2000;

phi_uvp_375 = zeros(sum(mask_uvp), 1);
for m = 1:numel(ia)
    id_uvp = ib(m);
    [~, iz] = min(abs(uvpd.depth_m(mask_z) - 375));
    phi_u = squeeze(uvpd.phi(id_uvp, mask_z, mask_uvp));
    if size(phi_u, 1) < size(phi_u, 2), phi_u = phi_u'; end
    phi_uvp_375 = phi_uvp_375 + phi_u(iz, :)';
end
phi_uvp_375 = phi_uvp_375 / numel(ia);

% ---------------------------------------------------------------
% 2. Three cases
% ---------------------------------------------------------------
cases = { ...
    struct('mode', 'off',     'label', 'DVM off',      'color', 'k'), ...
    struct('mode', 'gut',     'label', 'Gut only',     'color', 'r'), ...
    struct('mode', 'bianchi', 'label', 'Bianchi-lite', 'color', 'b')  ...
};

ratio_all = zeros(numel(k_plot), numel(cases));
spec_375  = zeros(numel(d_um_mod), numel(cases));
total_bv  = zeros(numel(cases), 1);

for ic = 1:numel(cases)
    cfg = cfg_base(cases{ic}.mode);
    cfg.validate();

    sim = ColumnSimulation(cfg, col_grid, prof);
    sim.rhs.G_gut = 0;
    sim.rhs.B_mig = 0;
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);

    % spinup
    for icyc = 1:max_cycles
        phi_before = mean(sum(Y + Yfp, 2));
        for i_day = 1:n_days
            sim.rhs.profile.eps = keps_day.eps(:, i_day);
            for i_step = 1:steps_per_day
                Y(k_bc, :) = phi_bc_daily(i_day, :);
                [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
                Y(k_bc, :) = phi_bc_daily(i_day, :);
            end
        end
        phi_after  = mean(sum(Y + Yfp, 2));
        rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
        if rel_change < spinup_tol
            fprintf('%s: converged at cycle %d\n', cases{ic}.label, icyc);
            break;
        end
    end

    % final run
    sim.rhs.G_gut = 0;
    sim.rhs.B_mig = 0;
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);
    phi_mod  = zeros(numel(k_plot), 1);
    spec_sum = zeros(cfg.n_sections, 1);
    n_cast   = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            Y(k_bc, :) = phi_bc_daily(i_day, :);
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
            Y(k_bc, :) = phi_bc_daily(i_day, :);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            Ytot    = Y + Yfp;
            phi_mod = phi_mod + sum(Ytot(k_plot, :), 2);
            [~, k375_mod] = min(abs(col_grid.z_centers - 375));
            spec_sum = spec_sum + Ytot(k375_mod, :)';
            n_cast   = n_cast + 1;
        end
    end

    phi_mod        = phi_mod / max(n_cast, 1);
    spec_375(:,ic) = spec_sum / max(n_cast, 1);
    ratio_all(:,ic)= phi_mod ./ phi_uvp_ref;
    total_bv(ic)   = sum(Y(:) + Yfp(:));
end

% ---------------------------------------------------------------
% 3. Print table
% ---------------------------------------------------------------
fprintf('\n--- BV ratio (model / UVP), 100-2000 um ---\n');
fprintf('%-8s  %-10s  %-10s  %-12s\n', 'Depth', 'DVM off', 'Gut only', 'Bianchi-lite');
for ki = 1:numel(k_plot)
    fprintf('%5.0f m    %5.2f      %5.2f      %5.2f\n', ...
        z_mod(ki), ratio_all(ki,1), ratio_all(ki,2), ratio_all(ki,3));
end

fprintf('\n--- Mass conservation / standing stock ---\n');
fprintf('DVM off:      total BV = %.4e  (0.00%% vs DVM off)\n', total_bv(1));
fprintf('Gut only:     total BV = %.4e  (%.2f%% vs DVM off)\n', total_bv(2), ...
    100 * (total_bv(2) - total_bv(1)) / max(total_bv(1), 1e-20));
fprintf('Bianchi-lite: total BV = %.4e  (%.2f%% vs DVM off)\n', total_bv(3), ...
    100 * (total_bv(3) - total_bv(1)) / max(total_bv(1), 1e-20));

% ---------------------------------------------------------------
% 4. Ratio plot
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 10 12], 'Color', 'white');
hold on;
for ic = 1:numel(cases)
    plot(ratio_all(:,ic), z_mod, [cases{ic}.color '-o'], ...
        'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end
plot([1 1], [z_mod(1) z_mod(end)], 'k:', 'LineWidth', 0.8, 'HandleVisibility', 'off');
set(gca, 'YDir', 'reverse', 'YLim', [60 510], 'XLim', [0 2.2]);
xlabel('Model / UVP ratio');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 7);
title('Bianchi-lite DVM test', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'bianchi_lite_ratio.png'));
fprintf('Saved bianchi_lite_ratio.png\n');

% ---------------------------------------------------------------
% 5. Spectrum at 375 m
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 11 8], 'Color', 'white');
hold on;
plot(uvpd.d_um(mask_uvp), phi_uvp_375, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', 'UVP');
for ic = 1:numel(cases)
    plot(d_um_mod(mask_mod), spec_375(mask_mod, ic), [cases{ic}.color '-s'], ...
        'MarkerSize', 3, 'LineWidth', 1.2, 'DisplayName', cases{ic}.label);
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('ESD (\mum)');
ylabel('\phi (ppmV)');
legend('Location', 'southwest', 'FontSize', 7);
title('375 m size spectrum', 'FontWeight', 'normal');
saveas(gcf, fullfile(fig_dir, 'bianchi_lite_spectrum_375m.png'));
fprintf('Saved bianchi_lite_spectrum_375m.png\n');

% ---------------------------------------------------------------
function cfg = cfg_base(mode)
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;   % dummy for validate; actual Dmax comes from Dmax_A
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_coag    = true;
cfg.enable_disagg  = true;
cfg.enable_zoo     = true;
cfg.enable_microbe = false;
cfg.enable_mining  = true;
cfg.alpha          = 0.10;
cfg.microbe_r0     = 0.0;
cfg.surface_pp_mu  = 0.0;
cfg.r_to_rg        = 1.6;
cfg.zoo_c          = 0.025;
cfg.zoo_s          = 1.3e-5;
cfg.zoo_p          = 0.5;
cfg.zoo_ic         = 7;
cfg.mining_s       = 1.3e-5;
cfg.fp_alpha_cross = 0.5;

cfg.enable_dvm          = false;
cfg.enable_dvm_bianchi  = false;
cfg.dvm_p               = 1.0;
cfg.dvm_ffec            = 0.2;
cfg.dvm_feed_zmax       = 150;
cfg.dvm_zmin            = 300;
cfg.dvm_zmax            = 500;
cfg.dvm_tau_gut         = 0.25;
cfg.dvm_tau_body        = 1.0;
cfg.dvm_body_frac       = 0.3;
cfg.dvm_mort_rate       = 0.02;
cfg.dvm_mort_bin        = 15;

if strcmpi(mode, 'gut')
    cfg.enable_dvm = true;
elseif strcmpi(mode, 'bianchi')
    cfg.enable_dvm_bianchi = true;
end
end

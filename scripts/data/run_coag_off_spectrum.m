% run_coag_off_spectrum.m
%
% Does the U-shape at 475 m come from coagulation?
%
% Run two cases (coag on / coag off) with identical flux BC at 100 m.
% Plot size spectrum at 475 m for both + UVP on the same axes.
%
% If the U-shape disappears with coag off -> coagulation is the cause.
% If it stays -> something else (sinking, disagg, etc.).
%
% Saves: docs/figures/coag_off_spectrum_475m.png

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- grid and spinup settings ---
col_grid      = ColumnGrid(1000, 20);
keps_day      = load_keps_daily(mat_path, col_grid.z_centers);
prof          = load_keps(mat_path, col_grid.z_centers);
z_centers     = col_grid.z_centers;
dz            = col_grid.dz;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 80;
k_bc          = 2;
k_plot        = 2:10;

% depth of interest: 475 m
z_target = 475;
[~, k_z] = min(abs(z_centers - z_target));

% --- load BC and UVP ---
cfg0         = cfg_coag(true);
bc           = get_daily_bc_at_depth(uvp_file, cfg0, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% sinking speed per model bin (for flux BC)
sim_tmp = ColumnSimulation(cfg0, col_grid, prof);
d_cm    = sim_tmp.size_grid.dcomb(:)';
d_um    = d_cm * 1e4;
w_bin   = 66 * d_cm .^ 0.62;

% --- UVP spectrum at 475 m, cast-day average ---
mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ~, ib] = intersect(bc.dates, uvpd.dates);

d_uvp_plot   = uvpd.d_um(mask_uvp);
phi_uvp_475  = zeros(1, sum(mask_uvp));

for m = 1:numel(ib)
    phi_u = squeeze(uvpd.phi(ib(m), mask_z, mask_uvp));
    if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
    [~, iz] = min(abs(uvpd.depth_m(mask_z) - z_target));
    phi_uvp_475 = phi_uvp_475 + phi_u(iz,:);
end
phi_uvp_475 = phi_uvp_475 / max(numel(ib), 1);

% --- run two cases ---
cases = {struct('coag', true,  'label', 'coag on',  'ls', 'k-'), ...
         struct('coag', false, 'label', 'coag off', 'ls', 'r--')};

phi_mod_475 = zeros(2, cfg0.n_sections);

for ic = 1:2
    cfg = cfg_coag(cases{ic}.coag);
    sim = ColumnSimulation(cfg, col_grid, prof);
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
            fprintf('%s: converged at cycle %d\n', cases{ic}.label, icyc);
            break;
        end
    end

    % accumulate spectrum at 475 m on cast days
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);
    n_cast = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc,:) = Y(k_bc,:) + flux_src;
            [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            phi_mod_475(ic,:) = phi_mod_475(ic,:) + (Y(k_z,:) + Yfp(k_z,:));
            n_cast = n_cast + 1;
        end
    end
    phi_mod_475(ic,:) = phi_mod_475(ic,:) / max(n_cast, 1);
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 8],'Color','white');
hold on;
plot(d_uvp_plot, phi_uvp_475, 'b-o', 'MarkerSize', 3, 'LineWidth', 1.2, ...
    'DisplayName', 'UVP');
for ic = 1:2
    plot(d_um, phi_mod_475(ic,:), cases{ic}.ls, 'LineWidth', 1.2, ...
        'DisplayName', cases{ic}.label);
end
set(gca,'XScale','log','YScale','log','XLim',[100 2000],'FontSize',fs,'Box','off');
xlabel('ESD (\mum)', 'FontSize', fs);
ylabel('BV (m^3 m^{-3})', 'FontSize', fs);
legend('Location','northeast','FontSize',fs,'Box','off');
title('475 m spectrum: coag on vs off', 'FontWeight','normal','FontSize',fs);

saveas(gcf, fullfile(fig_dir, 'coag_off_spectrum_475m.png'));
fprintf('Saved coag_off_spectrum_475m.png\n');

% ---------------------------------------------------------------
function cfg = cfg_coag(enable_coag)
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
cfg.disagg_dmax_A  = 9.39e-6 * 5;
cfg.enable_coag    = enable_coag;
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
end

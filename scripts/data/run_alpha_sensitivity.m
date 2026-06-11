% run_alpha_sensitivity.m
%
% Test: how does alpha (stickiness) affect the size spectrum?
%
% Higher alpha -> more coagulation -> faster growth to large sizes.
% Currently alpha = 0.5. Try 0.1, 0.3, 0.5, 0.7, 1.0.
% Compare spectrum at 75 m and 200 m to UVP.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

cfg_base = SimulationConfig();
cfg_base.n_sections     = 30;
cfg_base.sinking_law    = 'kriest_8';
cfg_base.disagg_mode    = 'operator_split';
cfg_base.disagg_dmax_cm = 1.0;
cfg_base.enable_coag    = true;
cfg_base.enable_disagg  = true;
cfg_base.enable_zoo     = true;
cfg_base.enable_microbe = true;
cfg_base.enable_mining  = true;
cfg_base.microbe_r0     = 0.03;
cfg_base.microbe_use_temp = true;
cfg_base.microbe_tref_C = 20;
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
d_model_mm = (2 * r_cm * 1e4)' / 1000;
log_d      = log(d_model_mm);
log_bnd    = [log_d(1)-(log_d(2)-log_d(1))/2, ...
              (log_d(1:end-1)+log_d(2:end))/2, ...
              log_d(end)+(log_d(end)-log_d(end-1))/2];
dw_model_mm = diff(exp(log_bnd));   % 1 x n_sec

% UVP
uvpd = parse_uvp_daily(uvp_file);
uvp_bin_mask = uvpd.d_um >= 100 & uvpd.d_um < 2000;
d_uvp_mm     = uvpd.d_um(uvp_bin_mask) / 1000;
dw_filt      = uvpd.dw(uvp_bin_mask);

[cast_dates_matched, ia, ib] = intersect(daily.dates, uvpd.dates);
[~, best]  = max(sum(daily.phi(ia, :), 2));
id_model   = ia(best);
id_uvp     = ib(best);
lbl = num2str(cast_dates_matched(best));
lbl = [lbl(5:6) '-' lbl(7:8)];

% alpha values to test
alphas = [0.1, 0.3, 0.5, 0.7, 1.0];
colors = {'b', 'c', 'k', 'r', 'm'};
check_z = [75, 200];

figure('Units', 'centimeters', 'Position', [2 2 18 8]);

for p = 1:2
    subplot(1, 2, p);
    hold on;

    [~, iz_mod] = min(abs(col_grid.z_centers - check_z(p)));
    [~, iz_uvp] = min(abs(uvpd.depth_m - check_z(p)));

    % UVP spectrum
    phi_u = squeeze(uvpd.phi(id_uvp, iz_uvp, uvp_bin_mask));
    S_u   = reshape(phi_u, 1, []) ./ dw_filt(:)' * 1e9;
    S_u(S_u <= 0) = NaN;
    ok_u = ~isnan(S_u);
    loglog(d_uvp_mm(ok_u), S_u(ok_u), 'g--', 'LineWidth', 2, 'DisplayName', 'UVP');

    for ka = 1:numel(alphas)
        cfg = cfg_base;
        cfg.alpha = alphas(ka);
        fprintf('z=%dm  alpha=%.1f ...\n', check_z(p), alphas(ka));

        cfg.validate();
        sim = ColumnSimulation(cfg, col_grid, prof);

        Y   = zeros(n_z, cfg.n_sections);
        Yfp = zeros(n_z, cfg.n_sections);
        for icyc = 1:max_cycles
            phi_before = mean(sum(Y + Yfp, 3), 2);
            for i_day = 1:n_days
                for i_step = 1:steps_per_day
                    Y(1,:) = daily.phi(i_day,:);
                    [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
                    Y(1,:) = daily.phi(i_day,:);
                end
            end
            phi_after  = mean(sum(Y + Yfp, 3), 2);
            rel_change = max(abs(phi_after - phi_before) ./ max(phi_before, 1e-20));
            if rel_change < spinup_tol, break; end
        end
        phi_snap = [];
        for i_day = 1:n_days
            for i_step = 1:steps_per_day
                Y(1,:) = daily.phi(i_day,:);
                [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
                Y(1,:) = daily.phi(i_day,:);
            end
            if i_day == id_model
                phi_snap = Y + Yfp;
            end
        end

        phi_m = reshape(phi_snap(iz_mod, :), 1, []);
        S_m   = phi_m ./ dw_model_mm * 1e6;
        dm_ok = reshape(d_model_mm(S_m > 0 & isfinite(S_m)), 1, []);
        Sm_ok = reshape(S_m(S_m > 0 & isfinite(S_m)), 1, []);
        if ~isempty(dm_ok)
            loglog(dm_ok, Sm_ok, colors{ka}, 'LineWidth', 1.2, ...
                'DisplayName', sprintf('a=%.1f', alphas(ka)));
        end
    end

    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlim([0.05 10]);
    ylim([1e-2 1e3]);
    xlabel('ESD (mm)');
    ylabel('ppmV mm^{-1}');
    legend('location', 'southwest', 'FontSize', 6);
    title(sprintf('z = %d m  date = %s', check_z(p), lbl));
end

saveas(gcf, fullfile(fig_dir, 'alpha_sensitivity.png'));
fprintf('Saved alpha_sensitivity.png\n');

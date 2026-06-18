% run_report_fig3_dvm.m
%
% Figure 3 for combined report.
% DVM test: compare Model/UVP ratio with DVM off vs DVM on.
%
% Uses the same flux BC as Fig 1 (w*phi_UVP/dz at layer 2).
% DVM routes a fraction of fecal production from surface zone
% down to 300-500 m. We test whether this moves enough material
% to close the deep deficit.
%
% Result: ratio profiles are nearly identical -- DVM has < 2% effect at 475 m.
%
% Saves: report_fig3_dvm.png

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir  = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- setup ---
col_grid  = ColumnGrid(1000, 20);
keps_day  = load_keps_daily(mat_path, col_grid.z_centers);
prof      = load_keps(mat_path, col_grid.z_centers);
z_centers = col_grid.z_centers;
n_z       = col_grid.n_z;
dz        = col_grid.dz;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 80;
k_plot        = 2:10;
z_mod         = z_centers(k_plot);
k_bc          = 2;

% --- load UVP ---
cfg0 = cfg_base(false);
bc   = get_daily_bc_at_depth(uvp_file, cfg0, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

sim_tmp  = ColumnSimulation(cfg0, col_grid, prof);
d_cm     = sim_tmp.size_grid.dcomb(:)';
w_bin    = (66 * d_cm .^ 0.62);
mask_uvp = bc.d_model_um >= 100 & bc.d_model_um < 2000;

mask_uvp_raw = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z_uvp   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ~, ib]   = intersect(bc.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ib)
    phi_u = squeeze(uvpd.phi(ib(m), mask_z_uvp, mask_uvp_raw));
    if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z_uvp) - z_mod(ki)));
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz,:));
    end
end
phi_uvp_ref = phi_uvp_ref / max(numel(ib), 1);
mask_ok     = phi_uvp_ref >= 0.01 * max(phi_uvp_ref);

% --- run DVM off and DVM on ---
dvm_flags = {false, true};
labels    = {'DVM off', 'DVM on'};
ratio_out = zeros(numel(k_plot), 2);

for ir = 1:2
    cfg = cfg_base(dvm_flags{ir});
    sim = ColumnSimulation(cfg, col_grid, prof);
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);

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
            fprintf('%s: converged at cycle %d\n', labels{ir}, icyc); break;
        end
    end

    % collect on cast days
    Y   = zeros(n_z, cfg.n_sections);
    Yfp = zeros(n_z, cfg.n_sections);
    phi_mod = zeros(numel(k_plot), 1);
    n_cast  = 0;

    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc,:) = Y(k_bc,:) + flux_src;
            [Y, Yfp]  = sim.rhs.stepY(Y, dt, Yfp);
        end
        if any(bc.dates(i_day) == uvpd.dates)
            Ytot    = Y(k_plot,:) + Yfp(k_plot,:);
            phi_mod = phi_mod + sum(Ytot(:, mask_uvp), 2);
            n_cast  = n_cast + 1;
        end
    end
    phi_mod = phi_mod / max(n_cast, 1);
    ratio_out(:, ir) = NaN;
    ratio_out(mask_ok, ir) = phi_mod(mask_ok) ./ phi_uvp_ref(mask_ok);
end

% print ratio table
fprintf('\nDepth    DVM off   DVM on\n');
for ki = 1:numel(k_plot)
    if mask_ok(ki)
        fprintf('%5.0f m   %.3f     %.3f\n', z_mod(ki), ratio_out(ki,1), ratio_out(ki,2));
    end
end

% --- figure ---
fs = 7;
figure('Units','centimeters','Position',[2 2 9 10],'Color','white');
hold on;

fill([0.5 2 2 0.5],[min(z_mod(mask_ok)) min(z_mod(mask_ok)) ...
    max(z_mod(mask_ok)) max(z_mod(mask_ok))], ...
    [0.92 0.92 0.92],'EdgeColor','none','HandleVisibility','off');

plot([1 1],[min(z_mod(mask_ok)) max(z_mod(mask_ok))],'k:','LineWidth',0.8, ...
    'HandleVisibility','off');

plot(ratio_out(:,1), z_mod, 'k-o', 'MarkerSize',4,'LineWidth',1.2,'DisplayName','DVM off');
plot(ratio_out(:,2), z_mod, 'b-o', 'MarkerSize',4,'LineWidth',1.2,'DisplayName','DVM on');

set(gca,'YDir','reverse','XScale','log','YLim',[60 510],'XLim',[0.05 5], ...
    'FontSize',fs,'Box','off');
xlabel('Model / UVP','FontSize',fs);
ylabel('Depth (m)','FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('DVM test','FontWeight','normal','FontSize',fs);

saveas(gcf, fullfile(fig_dir, 'report_fig3_dvm.png'));
fprintf('Saved report_fig3_dvm.png\n');

% --- config ---
function cfg = cfg_base(enable_dvm)
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.disagg_mode    = 'operator_split';
cfg.disagg_dmax_cm = 1.0;
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
cfg.enable_dvm     = enable_dvm;
cfg.dvm_p          = 1.0;
cfg.dvm_ffec       = 0.2;
cfg.dvm_feed_zmax  = 150;
cfg.dvm_zmin       = 300;
cfg.dvm_zmax       = 500;
end

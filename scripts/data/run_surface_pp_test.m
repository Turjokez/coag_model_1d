% run_surface_pp_test.m
%
% Compare two BC modes against UVP biovolume profile:
%   Mode A: Flux BC at 100 m (w * phi_UVP injected at layer 2)
%   Mode B: Surface PP at 0 m (mu*phi growth at layer 1)
%
% Both use same physics (alpha=0.10, Da*5, zoo on, mining on).
% Ratio = model BV / UVP BV at each depth (100-2000 um).
%
% Saves: surface_pp_uvp_compare.png

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
fig_dir  = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% --- shared setup ---
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
k_plot        = 2:10;   % layers 2-10 = 75-475 m
z_mod         = z_centers(k_plot);

% --- load UVP BC struct (for flux BC mode and UVP reference) ---
cfg0 = cfg_base_fluxbc();
bc   = get_daily_bc_at_depth(uvp_file, cfg0, col_grid, 100, k_plot);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;

% bin sizes from flux BC config
cfg_tmp  = cfg_base_fluxbc();
sim_tmp  = ColumnSimulation(cfg_tmp, col_grid, prof);
d_cm     = sim_tmp.size_grid.dcomb(:)';
w_bin    = (66 * d_cm .^ 0.62);
mask_uvp = bc.d_model_um >= 100 & bc.d_model_um < 2000;

% --- UVP reference biovolume at each depth (cast-day average) ---
mask_uvp_raw = uvpd.d_um >= 100 & uvpd.d_um < 2000;
mask_z_uvp   = uvpd.depth_m >= 50 & uvpd.depth_m <= 510;
[~, ia, ib]  = intersect(bc.dates, uvpd.dates);

phi_uvp_ref = zeros(numel(k_plot), 1);
for m = 1:numel(ia)
    phi_u = squeeze(uvpd.phi(ib(m), mask_z_uvp, mask_uvp_raw));
    if size(phi_u,1) < size(phi_u,2), phi_u = phi_u'; end
    for ki = 1:numel(k_plot)
        [~, iz] = min(abs(uvpd.depth_m(mask_z_uvp) - z_mod(ki)));
        phi_uvp_ref(ki) = phi_uvp_ref(ki) + sum(phi_u(iz,:));
    end
end
phi_uvp_ref = phi_uvp_ref / numel(ia);
mask_ok     = phi_uvp_ref >= 0.01 * max(phi_uvp_ref);

% ---------------------------------------------------------------
% Mode A: Flux BC at 100 m
% ---------------------------------------------------------------
fprintf('\n--- Mode A: Flux BC ---\n');
cfg_a = cfg_base_fluxbc();
sim_a = ColumnSimulation(cfg_a, col_grid, prof);
Y     = zeros(n_z, cfg_a.n_sections);
Yfp   = zeros(n_z, cfg_a.n_sections);
k_bc  = 2;

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim_a.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc,:) = Y(k_bc,:) + flux_src;
            [Y, Yfp]  = sim_a.rhs.stepY(Y, dt, Yfp);
        end
    end
    phi_after = mean(sum(Y + Yfp, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('Converged at cycle %d\n', icyc); break;
    end
end

phi_a = zeros(numel(k_plot), 1);
n_cast = 0;
Y = zeros(n_z, cfg_a.n_sections); Yfp = zeros(n_z, cfg_a.n_sections);
for i_day = 1:n_days
    sim_a.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim_a.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        Ytot_a = Y(k_plot,:) + Yfp(k_plot,:);
        phi_a = phi_a + sum(Ytot_a(:, mask_uvp), 2);
        n_cast = n_cast + 1;
    end
end
phi_a = phi_a / max(n_cast, 1);

% ---------------------------------------------------------------
% Mode B: Surface PP
% ---------------------------------------------------------------
fprintf('\n--- Mode B: Surface PP (mu=0.1) ---\n');
cfg_b = cfg_base_fluxbc();
cfg_b.enable_surface_pp = true;
cfg_b.surface_pp_mu     = 0.1;
cfg_b.surface_pp_bin    = 1;

sim_b = ColumnSimulation(cfg_b, col_grid, prof);
Y     = zeros(n_z, cfg_b.n_sections);
Yfp   = zeros(n_z, cfg_b.n_sections);

% seed layer 1 bin 1 so mu*phi can grow from non-zero
Y(1,1) = 1e-8;

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim_b.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            [Y, Yfp] = sim_b.rhs.stepY(Y, dt, Yfp);
        end
    end
    phi_after = mean(sum(Y + Yfp, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        fprintf('Converged at cycle %d\n', icyc); break;
    end
end

% scale surface PP output to match UVP at shallowest valid depth
Ytot_b = Y(k_plot,:) + Yfp(k_plot,:);
phi_b_raw = sum(Ytot_b(:, mask_uvp), 2);
i_ref     = find(mask_ok, 1, 'first');
scale_b   = phi_uvp_ref(i_ref) / max(phi_b_raw(i_ref), 1e-30);
phi_b     = phi_b_raw * scale_b;
fprintf('Surface PP scale factor (to match UVP at %.0f m): %.3f\n', ...
    z_mod(i_ref), scale_b);

% ---------------------------------------------------------------
% Martin b for both modes (BV flux, layers > 50 m)
% ---------------------------------------------------------------
k_100 = find(z_centers >= 100, 1);
k_500 = find(z_centers >= 500, 1);

% Mode A BV flux
Y = zeros(n_z, cfg_a.n_sections); Yfp = zeros(n_z, cfg_a.n_sections);
F_bv_a = zeros(n_z,1); n_cast = 0;
for i_day = 1:n_days
    sim_a.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day,:)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc,:) = Y(k_bc,:) + flux_src;
        [Y, Yfp]  = sim_a.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        for k = 1:n_z
            F_bv_a(k) = F_bv_a(k) + sum(w_bin .* Y(k,:));
        end
        n_cast = n_cast + 1;
    end
end
F_bv_a = F_bv_a / max(n_cast,1);

% Mode B BV flux (from converged Y)
F_bv_b = zeros(n_z,1);
for k = 1:n_z
    F_bv_b(k) = sum(w_bin .* Y(k,:)) * scale_b;
end
% recompute from converged sim_b state
Y_b = zeros(n_z, cfg_b.n_sections); Yfp_b = zeros(n_z, cfg_b.n_sections);
Y_b(1,1) = 1e-8;
for icyc = 1:max_cycles
    for i_day = 1:n_days
        sim_b.rhs.profile.eps = keps_day.eps(:, i_day);
        for i_step = 1:steps_per_day
            [Y_b, Yfp_b] = sim_b.rhs.stepY(Y_b, dt, Yfp_b);
        end
    end
end
for k = 1:n_z
    F_bv_b(k) = sum(w_bin .* Y_b(k,:)) * scale_b;
end

b_a = -log(F_bv_a(k_500)/F_bv_a(k_100)) / log(z_centers(k_500)/z_centers(k_100));
if F_bv_b(k_100) > 0 && F_bv_b(k_500) > 0
    b_b = -log(F_bv_b(k_500)/F_bv_b(k_100)) / log(z_centers(k_500)/z_centers(k_100));
else
    b_b = NaN;
end
fprintf('\nMartin b — Flux BC: %.3f   Surface PP: %.3f   Canonical: 0.858\n', b_a, b_b);

% ---------------------------------------------------------------
% Ratios
% ---------------------------------------------------------------
ratio_a = NaN(numel(k_plot),1);
ratio_b = NaN(numel(k_plot),1);
ratio_a(mask_ok) = phi_a(mask_ok) ./ phi_uvp_ref(mask_ok);
ratio_b(mask_ok) = phi_b(mask_ok) ./ phi_uvp_ref(mask_ok);

% ---------------------------------------------------------------
% Figure: 2 panels
% ---------------------------------------------------------------
figure('Units','centimeters','Position',[2 2 14 10],'Color','white');

subplot(1,2,1);
hold on;
plot(phi_uvp_ref, z_mod, 'b-o', 'MarkerSize',3,'LineWidth',1.2,'DisplayName','UVP');
plot(phi_a,       z_mod, 'k-',  'MarkerSize',3,'LineWidth',1.2,'DisplayName','Flux BC');
plot(phi_b,       z_mod, 'k--', 'MarkerSize',3,'LineWidth',1.2,'DisplayName','Surface PP');
set(gca,'YDir','reverse','XScale','log','YLim',[60 510],'FontSize',7);
xlabel('BV (m^3 m^{-3})');  ylabel('Depth (m)');
legend('Location','southeast','FontSize',6);
title('BV profile (100–2000 \mum)','FontWeight','normal');

subplot(1,2,2);
hold on;
plot(ratio_a, z_mod, 'k-o',  'MarkerSize',3,'LineWidth',1.2,'DisplayName', ...
    sprintf('Flux BC (b=%.2f)',b_a));
plot(ratio_b, z_mod, 'k--o', 'MarkerSize',3,'LineWidth',1.2,'DisplayName', ...
    sprintf('Surface PP (b=%.2f)',b_b));
plot([1 1],[60 510],'r:','LineWidth',1.0);
set(gca,'YDir','reverse','XScale','log','YLim',[60 510],'XLim',[0.05 10],'FontSize',7);
xlabel('Model / UVP');  ylabel('Depth (m)');
legend('Location','southeast','FontSize',6);
title('Ratio vs depth','FontWeight','normal');

saveas(gcf, fullfile(fig_dir, 'surface_pp_uvp_compare.png'));
fprintf('Saved surface_pp_uvp_compare.png\n');

% ---------------------------------------------------------------
function cfg = cfg_base_fluxbc()
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
end

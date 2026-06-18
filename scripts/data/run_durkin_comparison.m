% run_durkin_comparison.m
%
% Compare model size-resolved particle flux to Durkin gel trap data.
%
% Model flux:  F(z, d) = w(d) * Y(z, d) / V(d)   [particles m-2 d-1]
%   w(d)  = Kriest_8 sinking speed [m/day]
%   Y(z,d)= model biovolume concentration [m3 BV / m3 water]
%   V(d)  = pi/6 * d_m^3   particle volume [m3/particle]
%
% Durkin flux: particles m-2 d-1 from gel traps at 125, 330, 500 m
%   - flux_agg: aggregate particles (ID2)
%   - flux_fp:  fecal pellets (ID3 + ID6 + ID7)
%
% Size filter: 100-2000 um for both model and Durkin.
% Uses flux BC (more physical than Dirichlet for flux interpretation).

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');
trap_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'sediment_trap_durkin', 'raw', ...
    'cb6a494508_EXPORTS_EXPORTSNA_JC214_classified_geltrap_particlefluxes.sb');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% ---------------------------------------------------------------
% 1. Load Durkin trap data
% ---------------------------------------------------------------
trap = load_durkin_flux(trap_file);
fprintf('Trap depths: '); fprintf('%d ', trap.depths); fprintf('m\n');
fprintf('Trap bins:   '); fprintf('%d ', trap.d_um);   fprintf('um\n');

% filter Durkin to 100-2000 um
mask_trap = trap.d_um >= 100 & trap.d_um < 2000;
d_trap    = trap.d_um(mask_trap);
flux_agg  = trap.flux_agg(:, mask_trap);   % n_dep x n_bin_trap
flux_fp   = trap.flux_fp(:, mask_trap);

% comparison depths: pick 125, ~330, ~500 m
target_z = [125, 330, 500];
id_trap  = zeros(1, 3);
for i = 1:3
    [~, id_trap(i)] = min(abs(trap.depths - target_z(i)));
end
z_used = trap.depths(id_trap);
fprintf('\nUsing trap depths: %d %d %d m\n', z_used(1), z_used(2), z_used(3));

% ---------------------------------------------------------------
% 2. Run model (flux BC, best config)
% ---------------------------------------------------------------
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

k_bc   = 2;
dz     = col_grid.dz;
n_z    = col_grid.n_z;

dt            = 0.25;
steps_per_day = round(1 / dt);
spinup_tol    = 0.01;
max_cycles    = 50;

cfg = cfg_best();
k_plot_bc = 2:10;
bc  = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, k_plot_bc);
phi_bc_daily = bc.phi_bc_daily;
n_days       = bc.n_days;
uvpd         = bc.uvpd;
d_model_um   = bc.d_model_um;   % model bin centers [um]

% sinking speeds [m/day] and particle volumes [m3]
d_cm  = d_model_um * 1e-4;
d_m   = d_model_um * 1e-6;
w_bin = (66 * d_cm .^ 0.62)';     % 1 x n_sec [m/day]
V_bin = (pi/6) * d_m .^ 3;        % 1 x n_sec [m3/particle]
V_bin = V_bin';                    % 1 x n_sec

% layer indices for comparison depths
k_comp = zeros(1, 3);
for i = 1:3
    [~, k_comp(i)] = min(abs(col_grid.z_centers - z_used(i)));
end
fprintf('Model layers used: k=%d(%dm) k=%d(%dm) k=%d(%dm)\n', ...
    k_comp(1), col_grid.z_centers(k_comp(1)), ...
    k_comp(2), col_grid.z_centers(k_comp(2)), ...
    k_comp(3), col_grid.z_centers(k_comp(3)));

sim = ColumnSimulation(cfg, col_grid, prof);
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);

% spinup
for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, i_day);
        flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
        for i_step = 1:steps_per_day
            Y(k_bc, :) = Y(k_bc, :) + flux_src;
            [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
        end
    end
    phi_after  = mean(sum(Y + Yfp, 2));
    rel_change = abs(phi_after - phi_before) / max(phi_before, 1e-20);
    if rel_change < spinup_tol
        fprintf('Converged at cycle %d\n', icyc);
        break;
    end
end

% final run: accumulate model flux on cast days
Y   = zeros(n_z, cfg.n_sections);
Yfp = zeros(n_z, cfg.n_sections);
Y_acc   = zeros(3, cfg.n_sections);   % marine snow
Yfp_acc = zeros(3, cfg.n_sections);   % fecal pellets
n_cast  = 0;

for i_day = 1:n_days
    sim.rhs.profile.eps = keps_day.eps(:, i_day);
    flux_src = dt * (w_bin .* phi_bc_daily(i_day, :)) / dz;
    for i_step = 1:steps_per_day
        Y(k_bc, :) = Y(k_bc, :) + flux_src;
        [Y, Yfp]   = sim.rhs.stepY(Y, dt, Yfp);
    end
    if any(bc.dates(i_day) == uvpd.dates)
        for i = 1:3
            Y_acc(i, :)   = Y_acc(i, :)   + Y(k_comp(i), :);
            Yfp_acc(i, :) = Yfp_acc(i, :) + Yfp(k_comp(i), :);
        end
        n_cast = n_cast + 1;
    end
end

Y_mean   = Y_acc   / max(n_cast, 1);    % 3 x n_sec, mean BV conc
Yfp_mean = Yfp_acc / max(n_cast, 1);

% convert BV concentration -> number flux [particles m-2 d-1]
% F_n = w * Y / V_particle
mask_mod = d_model_um >= 100 & d_model_um < 2000;
F_agg = zeros(3, cfg.n_sections);    % marine snow number flux
F_fp  = zeros(3, cfg.n_sections);    % fecal pellet number flux
for i = 1:3
    F_agg(i, :) = w_bin .* Y_mean(i, :)   ./ V_bin;
    F_fp(i, :)  = w_bin .* Yfp_mean(i, :) ./ V_bin;
end

% ---------------------------------------------------------------
% 3. Print ratio table (total flux in valid size range)
% ---------------------------------------------------------------
fprintf('\n--- Total particle flux 100-2000 um [particles m-2 d-1] ---\n');
fprintf('%-10s  %-12s  %-12s  %-12s  %-8s\n', ...
    'Depth', 'Model agg', 'Model fp', 'Trap agg', 'Ratio');
for i = 1:3
    mod_tot  = sum(F_agg(i, mask_mod));
    modfp_tot = sum(F_fp(i, mask_mod));
    trap_tot = sum(flux_agg(id_trap(i), :), 'omitnan');
    ratio    = mod_tot / max(trap_tot, 1e-30);
    fprintf('%5.0f m    %10.2e    %10.2e    %10.2e    %6.2f\n', ...
        z_used(i), mod_tot, modfp_tot, trap_tot, ratio);
end

% ---------------------------------------------------------------
% 4. Plot: size spectrum at 3 depths
% ---------------------------------------------------------------
figure('Units', 'centimeters', 'Position', [2 2 18 6], 'Color', 'white');

for i = 1:3
    subplot(1, 3, i);
    hold on;

    % model aggregate flux
    plot(d_model_um(mask_mod), F_agg(i, mask_mod), 'k-o', ...
         'MarkerSize', 3, 'LineWidth', 1.0, 'DisplayName', 'Model (agg)');

    % model fecal flux
    plot(d_model_um(mask_mod), F_fp(i, mask_mod), 'k--s', ...
         'MarkerSize', 3, 'LineWidth', 0.8, 'DisplayName', 'Model (fp)');

    % Durkin aggregate
    plot(d_trap, flux_agg(id_trap(i), :), 'b-o', ...
         'MarkerSize', 4, 'LineWidth', 1.0, 'DisplayName', 'Trap (agg)');

    % Durkin fecal
    plot(d_trap, flux_fp(id_trap(i), :), 'b--s', ...
         'MarkerSize', 4, 'LineWidth', 0.8, 'DisplayName', 'Trap (fp)');

    set(gca, 'XScale', 'log', 'YScale', 'log', ...
        'XLim', [90 2100], 'FontSize', 7);
    xlabel('Diameter (\mum)');
    if i == 1
        ylabel('Flux (particles m^{-2} d^{-1})');
        legend('Location', 'northeast', 'FontSize', 5);
    end
    title(sprintf('%d m', z_used(i)), 'FontWeight', 'normal');
    hold off;
end

saveas(gcf, fullfile(fig_dir, 'durkin_comparison.png'));
fprintf('\nSaved durkin_comparison.png\n');

% ---------------------------------------------------------------
% 5. Plot: depth profile of total flux (100-2000 um)
% ---------------------------------------------------------------
mod_total_flux = sum(F_agg(:, mask_mod), 2);
trap_total_agg = sum(flux_agg(id_trap, :), 2, 'omitnan');

figure('Units', 'centimeters', 'Position', [2 2 8 12], 'Color', 'white');
hold on;
plot(mod_total_flux,  z_used, 'k-o', 'MarkerSize', 4, 'LineWidth', 1.2, ...
     'DisplayName', 'Model (agg)');
plot(trap_total_agg, z_used, 'b-s', 'MarkerSize', 5, 'LineWidth', 1.2, ...
     'DisplayName', 'Trap (agg)');
set(gca, 'YDir', 'reverse', 'XScale', 'log', ...
    'YLim', [100 550], 'FontSize', 7);
xlabel('Total flux (particles m^{-2} d^{-1})');
ylabel('Depth (m)');
legend('Location', 'southeast', 'FontSize', 7);
title('Total flux profile', 'FontWeight', 'normal');
hold off;

saveas(gcf, fullfile(fig_dir, 'durkin_flux_profile.png'));
fprintf('Saved durkin_flux_profile.png\n');

% ---------------------------------------------------------------
function cfg = cfg_best()
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

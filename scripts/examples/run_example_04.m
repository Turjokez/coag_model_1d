% run_example_04.m
% Example 4: microbial remineralization sensitivity.
%
% Compare flux depth profiles for r0 = 0, 0.01, 0.02 day^-1.
% All other settings match Example 3 (zoo + disagg + mining).
% Each run is scaled so the 100-m flux matches Th-234 = 4.65 mmol/m2/day.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));
addpath(fullfile(script_dir, '..', 'data'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

% --- paths ---
mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% --- grid and profiles ---
col_grid = ColumnGrid(1000, 20);
keps_day = load_keps_daily(mat_path, col_grid.z_centers);
prof     = load_keps(mat_path, col_grid.z_centers);

% --- base config (production setup, same as Example 3) ---
cfg_base = SimulationConfig();
cfg_base.n_sections     = 30;
cfg_base.sinking_law    = 'kriest_8';
cfg_base.ds_kernel_mode = 'sinking_law';
cfg_base.r_to_rg        = 1.6;
cfg_base.alpha          = 0.093;     % best-fit value
cfg_base.enable_coag    = true;
cfg_base.enable_disagg  = true;
cfg_base.disagg_mode    = 'operator_split';
cfg_base.disagg_dmax_A  = 9.39e-6 * 5;
cfg_base.enable_zoo     = true;
cfg_base.zoo_c          = 0.025;
cfg_base.zoo_s          = 1.3e-5;
cfg_base.zoo_p          = 0.5;
cfg_base.zoo_ic         = 7;
cfg_base.enable_mining  = true;

% --- BC at 100 m (bc_scale = 0.42, best fit) ---
bc           = get_daily_bc_at_depth(uvp_file, cfg_base, col_grid, 100, 2:10);
phi_bc_daily = bc.phi_bc_daily * 0.42;   % apply best-fit bc_scale
n_days       = bc.n_days;

% --- experiment loop ---
r0_vals  = [0,    0.01,  0.02];
labels   = {'r_0 = 0', 'r_0 = 0.01 d^{-1}', 'r_0 = 0.02 d^{-1}'};
ls_style = {'k-', 'k--', 'k:'};

F_th234 = 4.65;  % Th-234 flux at 100 m [mmol/m2/day] (reference scale)

z_c  = col_grid.z_centers;
dz   = col_grid.dz;
flux = zeros(col_grid.n_z, 3);

dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;

spinup_tol = 0.01;
max_cycles  = 80;

for ir = 1:3
    cfg = copy(cfg_base);
    cfg.enable_microbe = (r0_vals(ir) > 0);
    cfg.microbe_r0     = r0_vals(ir);

    sim   = ColumnSimulation(cfg, col_grid, prof);
    w_bin = 66 * sim.size_grid.dcomb(:)' .^ 0.62;

    % spinup
    Y   = zeros(col_grid.n_z, cfg.n_sections);
    Yfp = zeros(col_grid.n_z, cfg.n_sections);
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
        phi_after  = mean(sum(Y + Yfp, 2));
        if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
            fprintf('r0=%.3f  Converged cycle %d\n', r0_vals(ir), icyc);
            break;
        end
    end

    % total biovolume flux at each depth [m3/m2/day]
    for iz = 1:col_grid.n_z
        flux(iz, ir) = sum((Y(iz,:) + Yfp(iz,:)) .* w_bin);
    end

    % scale so 100-m flux matches Th-234
    k100 = find(z_c >= 100, 1, 'first');
    flux(:, ir) = flux(:, ir) * (F_th234 / flux(k100, ir));
end

% --- figure ---
fs = 8;
figure('Units','centimeters','Position',[2 2 8 9],'Color','white');
hold on;
for ir = 1:3
    plot(flux(:, ir), z_c, ls_style{ir}, 'LineWidth', 1.2, ...
         'DisplayName', labels{ir});
end
set(gca, 'YDir','reverse', 'FontSize',fs, 'Box','on', 'YLim',[0 1000]);
xlabel('flux (mmol m^{-2} d^{-1})', 'FontSize',fs);
ylabel('depth (m)', 'FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('microbial r_0 sensitivity', 'FontWeight','normal','FontSize',fs);

out_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(out_dir, 'dir'), mkdir(out_dir); end
print(gcf, fullfile(out_dir,'example_04_flux.png'), '-dpng', '-r150');
fprintf('Saved example_04_flux.png\n');

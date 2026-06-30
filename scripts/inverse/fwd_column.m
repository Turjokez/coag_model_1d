function phi_out = fwd_column(params, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base)
% fwd_column  Run the column model and return spectra at obs_depths.
%
% params     = [alpha, zoo_c_scale]  (row or col vector)
% obs_depths = target depths in m, e.g. [125, 325, 475]
% phi_out    = numel(obs_depths) x n_sections matrix

alpha_try      = params(1);
zoo_c_scale    = params(2);

% copy base config and set trial parameters
cfg = copy(cfg_base);
cfg.alpha   = alpha_try;
cfg.zoo_c   = cfg_base.zoo_c * zoo_c_scale;

sim   = ColumnSimulation(cfg, col_grid, prof);
d_cm  = sim.size_grid.dcomb(:)';
w_bin = 66 * d_cm .^ 0.62;           % kriest_8 sinking [m/day]

dz            = col_grid.dz;
dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;                    % flux BC at 100 m (layer 2)
spinup_tol    = 0.01;
max_cycles    = 80;

Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);

% spinup until quasi-steady state
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
        break
    end
end

% extract spectra at requested depths
n_dep    = numel(obs_depths);
phi_out  = zeros(n_dep, cfg.n_sections);
z_c      = col_grid.z_centers;
for id = 1:n_dep
    [~, k] = min(abs(z_c - obs_depths(id)));
    phi_out(id,:) = Y(k,:) + Yfp(k,:);
end
end

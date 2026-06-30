function phi_out = fwd_column_v3(params, obs_depths, col_grid, keps_day, prof, phi_bc_daily, n_days, cfg_base)
% fwd_column_v3  Forward model: [alpha, zoo_c_scale, r0].
%
% params = [alpha, zoo_c_scale, r0]
%   alpha        : stickiness
%   zoo_c_scale  : multiplier on cfg_base.zoo_c clearance rate
%   r0           : microbial remineralization rate [day^-1]
%
% BC flux is unchanged (bc_scale = 1). Only physics params vary.
%
% phi_out = numel(obs_depths) x n_sections

alpha_try   = params(1);
zoo_c_sc    = params(2);
r0          = params(3);

if alpha_try <= 0 || alpha_try > 2 || zoo_c_sc <= 0 || zoo_c_sc > 20 || r0 < 0 || r0 > 1
    phi_out = zeros(numel(obs_depths), cfg_base.n_sections);
    return
end

cfg = copy(cfg_base);
cfg.alpha          = alpha_try;
cfg.zoo_c          = cfg_base.zoo_c * zoo_c_sc;
cfg.enable_microbe = (r0 > 0);
cfg.microbe_r0     = r0;

sim   = ColumnSimulation(cfg, col_grid, prof);
w_bin = 66 * sim.size_grid.dcomb(:)' .^ 0.62;

dz            = col_grid.dz;
dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 2;
spinup_tol    = 0.01;
max_cycles    = 80;

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
    phi_after = mean(sum(Y + Yfp, 2));
    if abs(phi_after - phi_before) / max(phi_before, 1e-20) < spinup_tol
        break
    end
end

z_c     = col_grid.z_centers;
n_dep   = numel(obs_depths);
phi_out = zeros(n_dep, cfg.n_sections);
for id = 1:n_dep
    [~, k] = min(abs(z_c - obs_depths(id)));
    phi_out(id,:) = Y(k,:) + Yfp(k,:);
end
end

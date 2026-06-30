function phi_out = fwd_column_surface(params, obs_depths, col_grid, keps_day, prof, phi_surf_daily, n_days, cfg_base)
% FWD_COLUMN_SURFACE  Forward model with surface BC (k=1, z=25m).
%
% params = [alpha, bc_scale, r0]
%   alpha    : stickiness
%   bc_scale : scale on surface flux injection
%   r0       : microbial remineralization rate [day^-1]
%
% Surface BC: flux = w * phi_surf * bc_scale / dz injected at k=1
% phi_surf_daily from get_daily_surface_phi (depth <= 5m, with power-law fill)
%
% phi_out = numel(obs_depths) x n_sections

alpha_try = params(1);
bc_scale  = params(2);
r0        = params(3);

if alpha_try<=0 || alpha_try>2 || bc_scale<=0 || bc_scale>2 || r0<0 || r0>1
    phi_out = zeros(numel(obs_depths), cfg_base.n_sections);
    return
end

cfg = copy(cfg_base);
cfg.alpha          = alpha_try;
cfg.enable_microbe = true;
cfg.microbe_r0     = r0;

sim   = ColumnSimulation(cfg, col_grid, prof);
d_cm  = sim.size_grid.dcomb(:)';
w_bin = 66 * d_cm .^ 0.62;   % kriest_8 [m/day]

dz            = col_grid.dz;
dt            = 0.25;
steps_per_day = round(1/dt);
k_bc          = 1;            % inject at top layer
spinup_tol    = 0.01;
max_cycles    = 80;

Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);

for icyc = 1:max_cycles
    phi_before = mean(sum(Y + Yfp, 2));
    for i_day = 1:n_days
        sim.rhs.profile.eps = keps_day.eps(:, min(i_day, size(keps_day.eps,2)));
        flux_src = dt * (w_bin .* phi_surf_daily(i_day,:)) * bc_scale / dz;
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

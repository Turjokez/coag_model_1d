function obs = load_uvp_obs(uvp_file, obs_depths)
% load_uvp_obs  Load UVP observations at target depths.
%
% Filters to 100-2000 um (particle range; removes zooplankton).
% Returns total biovolume at each obs depth.
%
% obs.bv_total   [numel(obs_depths) x 1]  total BV [cm^3/cm^3]
% obs.depth      [numel(obs_depths) x 1]  actual depth used [m]

uvp = parse_uvp(uvp_file);

% filter to 100-2000 um
mask = uvp.d_um >= 100 & uvp.d_um < 2000;

n_dep    = numel(obs_depths);
bv_total = zeros(n_dep, 1);
dep_used = zeros(n_dep, 1);

for id = 1:n_dep
    [~, k] = min(abs(uvp.depth_m - obs_depths(id)));
    dep_used(id)  = uvp.depth_m(k);
    phi_row       = uvp.phi(k, :);
    bv_total(id)  = sum(phi_row(mask), 'omitnan');
end

obs.bv_total = bv_total;
obs.depth    = dep_used;
end

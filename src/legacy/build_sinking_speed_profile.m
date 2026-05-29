function out = build_sinking_speed_profile(base_speed_m_s, temp_c, rho_kg_m3)
% build_sinking_speed_profile
% Build depth-dependent sinking profile from relative viscosity scaling.

base_speed_m_s = base_speed_m_s(:)';
temp_c = temp_c(:);
rho_kg_m3 = rho_kg_m3(:);

% If salinity is not given in this helper path, keep a typical value.
sal_psu = 35.0 .* ones(size(temp_c));
nu = water_kinematic_viscosity_simple(temp_c, sal_psu, rho_kg_m3);

scale = nu(1) ./ max(nu, realmin);
scale(~isfinite(scale)) = 1.0;
scale(scale <= 0) = 1.0;

nz = numel(temp_c);
ns = numel(base_speed_m_s);
speed_profile_m_s = zeros(nz, ns);
for is = 1:ns
    speed_profile_m_s(:, is) = base_speed_m_s(is) .* scale;
end

out = struct();
out.scale = scale;
out.nu_m2_s = nu;
out.speed_profile_m_s = speed_profile_m_s;
out.speed_profile_m_day = speed_profile_m_s .* 86400.0;
end


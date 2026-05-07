function nu_m2_s = water_kinematic_viscosity_simple(temp_c, sal_psu, rho_kg_m3)
% water_kinematic_viscosity_simple
% Short note:
% 1. dynamic viscosity from simple temperature fit
% 2. small salinity correction
% 3. divide by density to get kinematic viscosity

Tk = temp_c + 273.15;
mu_w = 2.414e-5 .* 10 .^ (247.8 ./ (Tk - 140.0)); % Pa s
mu = mu_w .* (1.0 + 0.0015 .* (sal_psu - 35.0));

nu_m2_s = mu ./ max(rho_kg_m3, realmin);
end


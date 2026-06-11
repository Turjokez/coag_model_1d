function prof = load_keps(mat_path, z_model)
% LOAD_KEPS  Build a DepthProfile from keps_for_dave.mat.
%
% Usage:
%   prof = load_keps(mat_path, z_model)
%
% Inputs:
%   mat_path  - full path to keps_for_dave.mat
%   z_model   - model depth cell centers [m], positive down, n_z x 1
%
% Output:
%   prof      - DepthProfile object with real T, S, rho, nu, eps, Kz
%
% Notes:
%   - z in mat file is negative (ocean convention). Flip to positive.
%   - eps units in mat file: m^2/s^3. Model uses cm^2/s^3. Multiply by 1e4.
%   - T, S, rho are time-mean over all valid (non-NaN) profiles.
%   - Kz is set from kappa_T (thermal diffusivity) time-mean.
%   - Profiles are interpolated to z_model using pchip, then extrapolated
%     as constant beyond the data range.

% load
raw = load(mat_path);
S   = raw.S;

% depth: negative -> positive, sort ascending
z_data = abs(S.z(:));              % m, positive, 300 x 1
[z_data, idx] = sort(z_data);     % sort shallow to deep

% time-mean each variable, ignoring missing values
eps_data  = mean_no_nan(S.eps(idx,:),  2);   % m^2/s^3
T_data    = mean_no_nan(S.T(idx,:),    2);   % deg C
Sal_data  = mean_no_nan(S.S(idx,:),    2);   % psu
rho_data  = mean_no_nan(S.rho(idx,:),  2);   % kg/m^3
kT_data   = mean_no_nan(S.kappa_T(idx,:), 2); % m^2/s

% unit conversions
% rho: kg/m^3 -> g/cm^3
rho_data  = rho_data * 1e-3;

% eps: m^2/s^3 -> cm^2/s^3
eps_data  = eps_data * 1e4;

% floor: set minimum eps to realistic deep-ocean background
% open-ocean below mixed layer: ~1e-8 W/kg = 1e-4 cm^2/s^3
% avoids the VMP noise floor (1e-12 m^2/s^3) giving unrealistically huge D_max
eps_floor = 1e-4;   % cm^2/s^3
eps_data  = max(eps_data, eps_floor);

% kinematic viscosity from temperature (Sharqawy 2010 approximation)
% nu [cm^2/s] ~ 0.01 * exp(0.02 * (20 - T_C))
nu_data = 0.01 * exp(0.02 * (20 - T_data));  % cm^2/s

% interpolate to model grid, clamp extrapolation to edge values
T_C_m   = interp_clamped(z_data, T_data,   z_model);
T_K_m   = T_C_m + 273.15;
Sal_m   = interp_clamped(z_data, Sal_data,  z_model);
rho_m   = interp_clamped(z_data, rho_data,  z_model);
nu_m    = interp_clamped(z_data, nu_data,   z_model);
eps_m   = interp_clamped(z_data, eps_data,  z_model);
Kz_m    = interp_clamped(z_data, kT_data,   z_model);

% Kz floor (at least 1e-5 m^2/s)
Kz_m = max(Kz_m, 1e-5);

% build DepthProfile
prof = DepthProfile(z_model(:), T_K_m(:), Sal_m(:), rho_m(:), ...
                    nu_m(:), eps_m(:), Kz_m(:));

% add Stemmann zoo profiles (same as typical())
[prof.Zc, prof.Zf, prof.Zm] = DepthProfile.stemmannZoo(z_model(:));

end

function y = mean_no_nan(x, dim)
% Mean without needing the Statistics Toolbox.
good = ~isnan(x);
x(~good) = 0;
n = sum(good, dim);
y = sum(x, dim) ./ max(n, 1);
y(n == 0) = NaN;
end

function yi = interp_clamped(z, f, zi)
% Interpolate and keep values inside the data range.
good = isfinite(z) & isfinite(f);
z = z(good);
f = f(good);
yi = interp1(z, f, zi, 'pchip', 'extrap');
yi = max(min(f), min(max(f), yi));
end

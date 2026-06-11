function kd = load_keps_daily(mat_path, z_model)
% LOAD_KEPS_DAILY  Daily-mean eps(z) profiles from keps_for_dave.mat.
%
% Usage:
%   kd = load_keps_daily(mat_path, z_model)
%
% Output:
%   kd.dates      - YYYYMMDD for each day [n_days x 1]
%   kd.eps        - daily-mean eps on model grid [n_z x n_days], cm^2/s^3
%   kd.eps_mean   - cruise-mean eps [n_z x 1] (same as load_keps result)
%
% Notes:
%   - eps below data range (< 300 m) held constant at deepest value.
%   - eps floor: 1e-8 cm^2/s^3 (avoids infinite D_max).

raw = load(mat_path);
S   = raw.S;

% depth: negative -> positive, sort shallow to deep
z_data = abs(S.z(:));
[z_data, idx] = sort(z_data);

% eps: 300 x 1293 in m^2/s^3 -> cm^2/s^3
eps_all = S.eps(idx, :) * 1e4;   % 300 x n_times

% mtime: MATLAB datenum -> YYYYMMDD
mtime = S.mtime(1, :);   % 1 x n_times
dates_all = datenum_to_yyyymmdd(mtime);   % 1 x n_times

% unique days in order
u_dates = unique(dates_all);
n_days  = numel(u_dates);
n_z     = numel(z_model);

eps_floor = 1e-4;   % cm^2/s^3  (realistic deep-ocean background)

eps_daily = zeros(n_z, n_days);
for d = 1:n_days
    mask = (dates_all == u_dates(d));
    % mean ignoring NaN
    e_day = mean_no_nan(eps_all(:, mask), 2);   % 300 x 1
    e_day = max(e_day, eps_floor);
    eps_daily(:, d) = interp_clamped(z_data, e_day, z_model);
end

kd.dates    = u_dates(:);
kd.eps      = eps_daily;
kd.eps_mean = mean(eps_daily, 2);

end

% -------------------------------------------------------------------------
function y = mean_no_nan(x, dim)
good = ~isnan(x);
x(~good) = 0;
n = sum(good, dim);
y = sum(x, dim) ./ max(n, 1);
y(n == 0) = NaN;
end

function yi = interp_clamped(z, f, zi)
good = isfinite(z) & isfinite(f);
yi = interp1(z(good), f(good), zi(:), 'pchip', 'extrap');
lo = min(f(good));  hi = max(f(good));
yi = max(lo, min(hi, yi));
end

function d = datenum_to_yyyymmdd(ml)
% Convert MATLAB datenum to YYYYMMDD integer array.
% MATLAB datenum 1 = Jan 1, year 1 (proleptic). Use datevec.
% Works via mod arithmetic to avoid toolbox dependency.
d = zeros(size(ml));
for i = 1:numel(ml)
    v = datevec(ml(i));
    d(i) = v(1)*10000 + v(2)*100 + v(3);
end
end

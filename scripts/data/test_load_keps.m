% test_load_keps.m
% Check that load_keps runs and produces sensible profiles.
%
% Steps:
%   1. Load keps_for_dave.mat and build DepthProfile
%   2. Check eps, T, nu are in valid ranges
%   3. Plot eps(z) vs synthetic profile from DepthProfile.typical()
%   4. Print PASS/FAIL for each check

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', '..', 'data', 'NA', 'Turbulance', 'keps_for_dave.mat');

% model grid: 1000 m, 20 layers (dz = 50 m)
z_model = (25 : 50 : 975)';   % 20 cell centers [m]

fprintf('=== test_load_keps ===\n');

% run loader
prof = load_keps(mat_path, z_model);

% compare to synthetic
prof_syn = DepthProfile.typical(z_model);

% ---- checks ----
pass = true;

% 1. eps in physical range [1e-9, 1e-1] cm^2/s^3
if all(prof.eps >= 1e-9) && all(prof.eps <= 1e-1)
    fprintf('PASS  eps range: %.2e to %.2e cm^2/s^3\n', min(prof.eps), max(prof.eps));
else
    fprintf('FAIL  eps out of range: %.2e to %.2e\n', min(prof.eps), max(prof.eps));
    pass = false;
end

% 2. T_K reasonable: 275 to 300 K (2 to 27 C)
T_C = prof.T_K - 273.15;
if all(T_C >= 2) && all(T_C <= 27)
    fprintf('PASS  T range: %.1f to %.1f C\n', min(T_C), max(T_C));
else
    fprintf('FAIL  T out of range: %.1f to %.1f C\n', min(T_C), max(T_C));
    pass = false;
end

% 3. nu in range [0.005, 0.02] cm^2/s
if all(prof.nu >= 0.005) && all(prof.nu <= 0.02)
    fprintf('PASS  nu range: %.4f to %.4f cm^2/s\n', min(prof.nu), max(prof.nu));
else
    fprintf('FAIL  nu out of range: %.4f to %.4f\n', min(prof.nu), max(prof.nu));
    pass = false;
end

% 4. eps decreases with depth on average (real data should show this trend)
eps_upper = mean(prof.eps(z_model < 150));
eps_lower = mean(prof.eps(z_model > 150));
if eps_upper > eps_lower
    fprintf('PASS  eps decreases with depth (upper %.2e > lower %.2e)\n', ...
        eps_upper, eps_lower);
else
    fprintf('WARN  eps does not decrease with depth - check data\n');
end

% 5. Zc, Zf, Zm set
if ~isempty(prof.Zc) && ~isempty(prof.Zf) && ~isempty(prof.Zm)
    fprintf('PASS  zoo profiles set (Zc, Zf, Zm)\n');
else
    fprintf('FAIL  zoo profiles missing\n');
    pass = false;
end

% ---- summary ----
if pass
    fprintf('ALL PASS\n');
else
    fprintf('SOME CHECKS FAILED\n');
end

% ---- plot: real vs synthetic eps(z) ----
figure;
semilogx(prof.eps,     z_model, 'b-o', 'MarkerSize', 4, 'DisplayName', 'keps data (mean)');
hold on;
semilogx(prof_syn.eps, z_model, 'r--',               'DisplayName', 'synthetic (typical)');
set(gca, 'YDir', 'reverse');
xlabel('\epsilon  [cm^2 s^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('eps(z): data vs synthetic');

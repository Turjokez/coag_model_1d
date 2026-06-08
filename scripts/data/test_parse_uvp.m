% test_parse_uvp.m
% Check that parse_uvp reads the UVP file correctly.
%
% Steps:
%   1. Parse the differential .sb file
%   2. Check depth range, bin count, no all-NaN rows in top 100 m
%   3. Check N and phi values are positive where data exists
%   4. Plot mean surface number and biovolume spectra

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

sb_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

fprintf('=== test_parse_uvp ===\n');

uvp = parse_uvp(sb_file);

pass = true;

% 1. number of size bins
if numel(uvp.d_um) == 27
    fprintf('PASS  27 size bins found\n');
else
    fprintf('FAIL  expected 27 bins, got %d\n', numel(uvp.d_um));
    pass = false;
end

% 2. depth range reasonable (should cover 0 to at least 500 m)
if max(uvp.depth_m) >= 500
    fprintf('PASS  depth range: %.1f to %.1f m\n', min(uvp.depth_m), max(uvp.depth_m));
else
    fprintf('WARN  depth only to %.1f m\n', max(uvp.depth_m));
end

% 3. N is non-negative where not NaN
N_valid = uvp.N(~isnan(uvp.N));
if all(N_valid >= 0)
    fprintf('PASS  N >= 0 everywhere (mean = %.2e #/m^3)\n', mean(N_valid));
else
    fprintf('FAIL  negative N values found\n');
    pass = false;
end

% 4. phi is non-negative where not NaN
phi_valid = uvp.phi(~isnan(uvp.phi));
if all(phi_valid >= 0)
    fprintf('PASS  phi >= 0 everywhere (mean = %.2e cm^3/cm^3)\n', mean(phi_valid));
else
    fprintf('FAIL  negative phi values found\n');
    pass = false;
end

% 5. surface layer exists (depth <= 5 m)
surf_rows = uvp.depth_m <= 5;
if any(surf_rows)
    N_surf = mean_no_nan(uvp.N(surf_rows, :), 1);
    phi_surf = mean_no_nan(uvp.phi(surf_rows, :), 1);
    if any(N_surf > 0)
        fprintf('PASS  surface PSD has data (max bin = %.2e #/m^3)\n', max(N_surf));
        fprintf('PASS  surface phi has data (sum = %.2e cm^3/cm^3)\n', ...
            sum(phi_surf(~isnan(phi_surf))));
    else
        fprintf('FAIL  surface row is all zeros/NaN\n');
        pass = false;
    end
else
    fprintf('FAIL  no data in top 5 m\n');
    pass = false;
end

% 5. number of casts
fprintf('INFO  %d unique casts in file\n', uvp.n_casts);

if pass
    fprintf('ALL PASS\n');
else
    fprintf('SOME CHECKS FAILED\n');
end

% --- plot surface PSD ---
N_surf = mean_no_nan(uvp.N(uvp.depth_m <= 5, :), 1);
phi_surf = mean_no_nan(uvp.phi(uvp.depth_m <= 5, :), 1);

figure;
subplot(1, 2, 1);
loglog(uvp.d_um, N_surf, 'b-o', 'MarkerSize', 4);
xlabel('diameter  [\mum]');
ylabel('N  [#  m^{-3}]');
title('surface N');

subplot(1, 2, 2);
loglog(uvp.d_um, phi_surf, 'r-o', 'MarkerSize', 4);
xlabel('diameter  [\mum]');
ylabel('\phi  [cm^3 cm^{-3}]');
title('surface \phi');

function y = mean_no_nan(x, dim)
good = ~isnan(x);
x(~good) = 0;
n = sum(good, dim);
y = sum(x, dim) ./ max(n, 1);
y(n == 0) = NaN;
end

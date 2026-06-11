% plot_dmax_profile.m
%
% Quick diagnostic: what D_max does each disagg formula give at each depth?
% No simulation needed — just compute from the eps(z) profile.
%
% Shows:
%   Left  : eps(z) from VMP (cruise mean)
%   Right : D_max(z) for hard cutoff and logistic formulas

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

mat_path = fullfile(script_dir, '..', '..', 'data', 'NA', ...
    'Turbulance', 'keps_for_dave.mat');
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');

col_grid = ColumnGrid(1000, 20);
prof     = load_keps(mat_path, col_grid.z_centers);

z    = col_grid.z_centers;   % depth [m], 20 layers
eps  = prof.eps;             % cm^2/s^3

% --- hard cutoff: D_max = Dmax_A * eps_m^(-1/4) ---
Dmax_A   = 9.39e-6;                  % m
eps_m    = eps / 1e4;                % m^2/s^3
dmax_hard_mm = 1000 * Dmax_A * eps_m.^(-1/4);   % mm

% --- logistic: r_max = C0 * eps^(-B), D_max = 2*r_max (in mm) ---
C0 = 2e-3;   % cm
B  = 0.45;
r_max_cm     = C0 * eps.^(-B);
dmax_logistic_mm = 2 * r_max_cm * 10;   % cm -> mm

% --- print table ---
fprintf('%8s  %10s  %12s  %14s\n', 'z (m)', 'eps (cm2/s3)', 'D_max hard (mm)', 'D_max logistic (mm)');
for k = 1:numel(z)
    fprintf('%8.0f  %10.2e  %12.3f  %14.3f\n', z(k), eps(k), dmax_hard_mm(k), dmax_logistic_mm(k));
end

% --- plot ---
figure('Units', 'centimeters', 'Position', [2 2 16 10]);

subplot(1, 2, 1);
semilogx(eps, z, 'k', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse');
xlabel('\epsilon (cm^2 s^{-3})');
ylabel('depth (m)');
title('turbulent dissipation');

subplot(1, 2, 2);
semilogx(dmax_hard_mm, z, 'b', 'LineWidth', 1.5, 'DisplayName', 'hard D\_max');
hold on;
semilogx(dmax_logistic_mm, z, 'r', 'LineWidth', 1.5, 'DisplayName', 'logistic');
xline(0.3, 'g--', 'UVP min', 'LabelVerticalAlignment','bottom', 'FontSize', 7);
xline(2.0, 'g:', 'UVP max', 'LabelVerticalAlignment','bottom', 'FontSize', 7);
set(gca, 'YDir', 'reverse');
xlabel('D_{max} (mm)');
ylabel('depth (m)');
legend('location', 'southeast', 'FontSize', 7);
title('max stable size');

saveas(gcf, fullfile(fig_dir, 'dmax_profile.png'));
fprintf('\nSaved dmax_profile.png\n');

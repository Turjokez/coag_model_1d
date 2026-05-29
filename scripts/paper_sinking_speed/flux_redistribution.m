% flux_redistribution.m
% Keep total mass fixed, change PSD slope, and see flux/speed change.

clear; close all; clc;

%% --- Parameters ---

d = logspace(log10(10e-6), log10(5e-3), 1000);  % diameter [m], 10 um to 5 mm

xi_range = linspace(2.5, 5.0, 60);   % slopes to test (typical ocean range)

% Sinking speed: w = a * d^b [m/day]
% b from Kriest-type marine snow scaling, w_ref picked at 1 mm.
b = 0.62;
d_ref = 1e-3;        % 1 mm reference size
w_ref = 50;          % m/day at d_ref
a = w_ref / (d_ref^b);

M_target   = 10e-6;   % 10 mg/m^3 = 1e-5 kg/m^3
rho_eff    = 1050;    % effective density [kg/m^3]

%% --- Loop over slopes ---

flux   = zeros(size(xi_range));
w_mean = zeros(size(xi_range));

for k = 1:length(xi_range)
    xi = xi_range(k);

    % unnormalized number distribution: N ~ d^(-xi)
    N = d .^ (-xi);

    % mass per unit size
    m = (pi/6) .* rho_eff .* d.^3 .* N;

    % normalize so total mass = M_target
    N0   = M_target / trapz(d, m);
    m    = m * N0;

    % sinking speed at each size
    w = a .* d .^ b;

    % total flux [kg m^-2 day^-1] and mean speed [m/day]
    flux(k)   = trapz(d, w .* m);
    w_mean(k) = flux(k) / M_target;
end

% convert to mg m^-2 day^-1
flux_plot = flux * 1e6;

%% --- Plot ---

t = tiledlayout(1,2);
nexttile;
plot(xi_range, flux_plot, 'k', 'LineWidth', 1.5)
xlabel('\xi (size distribution slope)')
ylabel('flux (mg m^{-2} d^{-1})')
title('Total flux')

nexttile;
plot(xi_range, w_mean, 'k', 'LineWidth', 1.5)
xlabel('\xi (size distribution slope)')
ylabel('mean sinking speed (m d^{-1})')
title('Mean sinking speed')
set(gcf, 'Color', 'w');

fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'flux_redistribution.png'));

% quick summary
fprintf('flux range:   %.1f -- %.1f mg/m^2/day\n', min(flux_plot), max(flux_plot));
fprintf('speed range:  %.1f -- %.1f m/day\n', min(w_mean), max(w_mean));
fprintf('flux factor:  %.2f x\n', max(flux_plot) / max(min(flux_plot), eps));
fprintf('speed factor: %.2f x\n', max(w_mean) / max(min(w_mean), eps));
fprintf('saved: docs/figures/flux_redistribution.png\n');

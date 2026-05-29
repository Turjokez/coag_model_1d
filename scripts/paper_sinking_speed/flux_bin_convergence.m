% flux_bin_convergence.m
% Check how many size bins are enough to reproduce continuous flux.

clear; close all; clc;

%% Parameters
d_min = 10e-6;   % 10 um [m]
d_max = 5e-3;    % 5 mm [m]
xi_range = linspace(2.5, 5.0, 60);
n_bins_list = [1 2 5 10 20 30];

% Sinking speed: w = a * d^b [m/day]
b = 0.62;
d_ref = 1e-3;
w_ref = 50;
a = w_ref / (d_ref^b);

M_target = 10e-6;  % 10 mg/m^3 [kg/m^3]
rho_eff  = 1050;   % [kg/m^3]

%% Continuous reference
d_ref_grid = logspace(log10(d_min), log10(d_max), 3000);
flux_cont = zeros(size(xi_range));
for k = 1:numel(xi_range)
    xi = xi_range(k);
    flux_cont(k) = compute_flux_continuous(d_ref_grid, xi, a, b, M_target, rho_eff);
end

%% Discrete-bin flux and error
flux_disc = zeros(numel(n_bins_list), numel(xi_range));
err_pct   = zeros(size(flux_disc));

for i = 1:numel(n_bins_list)
    nb = n_bins_list(i);
    edges = logspace(log10(d_min), log10(d_max), nb + 1);
    dmid  = sqrt(edges(1:end-1) .* edges(2:end));

    for k = 1:numel(xi_range)
        xi = xi_range(k);

        % constant from mass normalization
        int_m = int_pow(edges(end), 4 - xi) - int_pow(edges(1), 4 - xi);
        Cmass = (pi/6) * rho_eff;
        N0 = M_target / max(Cmass * int_m, eps);

        % bin mass and bin flux
        m_bin = zeros(nb, 1);
        for j = 1:nb
            int_bin = int_pow(edges(j+1), 4 - xi) - int_pow(edges(j), 4 - xi);
            m_bin(j) = Cmass * N0 * int_bin;
        end

        w_mid = a .* dmid .^ b;
        flux_disc(i, k) = sum(w_mid(:) .* m_bin(:));

        err_pct(i, k) = 100 * abs(flux_disc(i, k) - flux_cont(k)) / max(abs(flux_cont(k)), eps);
    end
end

%% Summary table
fprintf('\n%-8s %-14s %-14s\n', 'n_bins', 'mean_err(%)', 'max_err(%)');
for i = 1:numel(n_bins_list)
    fprintf('%-8d %-14.3f %-14.3f\n', n_bins_list(i), mean(err_pct(i,:)), max(err_pct(i,:)));
end

% smallest n with max error <= 20%
ok = find(max(err_pct, [], 2) <= 20, 1, 'first');
if ~isempty(ok)
    fprintf('\nFirst n_bins with max error <= 20%%: %d\n', n_bins_list(ok));
else
    fprintf('\nNo case reached max error <= 20%% in tested n_bins.\n');
end

%% Plots
figure;
tiledlayout(1,2);

nexttile;
plot(xi_range, flux_cont * 1e6, 'k-', 'LineWidth', 1.8); hold on;
clr = lines(numel(n_bins_list));
for i = 1:numel(n_bins_list)
    plot(xi_range, flux_disc(i,:) * 1e6, '-', 'Color', clr(i,:), 'LineWidth', 1.1);
end
hold off;
xlabel('\xi (size distribution slope)');
ylabel('flux (mg m^{-2} d^{-1})');
lgd = cell(1, numel(n_bins_list) + 1);
lgd{1} = 'continuous';
for i = 1:numel(n_bins_list), lgd{i+1} = sprintf('%d bins', n_bins_list(i)); end
legend(lgd, 'Location', 'best');
title('Continuous vs discrete flux');

nexttile;
for i = 1:numel(n_bins_list)
    plot(xi_range, err_pct(i,:), '-', 'Color', clr(i,:), 'LineWidth', 1.2); hold on;
end
hold off;
yline(20, 'k--');
xlabel('\xi (size distribution slope)');
ylabel('flux error (%)');
legend(arrayfun(@(x) sprintf('%d bins', x), n_bins_list, 'UniformOutput', false), 'Location', 'best');
title('Discrete error vs slope');

set(gcf, 'Color', 'w');

fig_dir = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'flux_bin_convergence.png'));
fprintf('saved: docs/figures/flux_bin_convergence.png\n');

%% Local helpers
function F = compute_flux_continuous(d, xi, a, b, M_target, rho_eff)
    N = d .^ (-xi);
    m = (pi/6) .* rho_eff .* d.^3 .* N;
    N0 = M_target / max(trapz(d, m), eps);
    m = m * N0;
    w = a .* d .^ b;
    F = trapz(d, w .* m);
end

function val = int_pow(x, p)
    % integral helper: integral of d^(p-1) from 0 to x equals x^p/p
    if abs(p) < 1e-12
        val = log(x);
    else
        val = x.^p / p;
    end
end

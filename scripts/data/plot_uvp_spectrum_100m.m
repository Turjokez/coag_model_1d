% plot_uvp_spectrum_100m.m
%
% Show UVP volume spectrum at 100 m with power-law extension to sub-100 um.
% Demonstrates the same kind of merged-instrument spectrum plots Adrian showed
% (Jackson & Burd 2000 Fig 5; Zhang et al. 2022 Fig 13) but using only UVP
% data and our power-law extrapolation for the model BC.
%
% y-axis: volume spectral density  phi / Delta_d  [ppmV um^-1]
% x-axis: diameter [um], log scale

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

% minimal config just to get model bin sizes
cfg = SimulationConfig();
cfg.n_sections  = 30;
cfg.sinking_law = 'kriest_8';
cfg.validate();
col_grid = ColumnGrid(1000, 20);

bc = get_daily_bc_at_depth(uvp_file, cfg, col_grid, 100, 3:10);

% ---------------------------------------------------------------
% 1. UVP spectrum at BC depth on best cast day
% ---------------------------------------------------------------
uvpd   = bc.uvpd;
d_uvp  = bc.d_uvp_ok;     % UVP bin centers [um], 100-2000 um
dw_uvp = bc.dw_uvp_ok;    % UVP bin widths  [um]
d_uvp  = d_uvp(:);
dw_uvp = dw_uvp(:);

mask_uvp = uvpd.d_um >= 100 & uvpd.d_um < 2000;
phi_uvp  = squeeze(uvpd.phi(bc.id_uvp_best, bc.iz_bc, mask_uvp));
phi_uvp(isnan(phi_uvp)) = 0;

spec_uvp = phi_uvp(:) ./ dw_uvp(:);   % volume spectral density

% ---------------------------------------------------------------
% 2. Power-law fit on 100-400 um (same range as BC fill)
% ---------------------------------------------------------------
fit_ok = d_uvp >= 100 & d_uvp <= 400 & spec_uvp > 0;
p      = polyfit(log10(d_uvp(fit_ok)), log10(spec_uvp(fit_ok)), 1);
slope  = p(1);
fprintf('Volume spectrum slope: %.2f  (Junge xi = %.2f)\n', slope, 3 - slope);

% extended line from 1 um to 3000 um
d_line    = logspace(0, log10(3000), 400);
spec_line = 10 .^ polyval(p, log10(d_line));

% ---------------------------------------------------------------
% 3. Model bin spectral densities from BC
% ---------------------------------------------------------------
d_model = bc.d_model_um;
n_sec   = cfg.n_sections;

% model bin edges (same as get_daily_bc_at_depth)
d_edges = zeros(1, n_sec + 1);
d_edges(1)        = d_model(1)^2  / d_model(2);
d_edges(n_sec+1)  = d_model(n_sec)^2 / d_model(n_sec-1);
for k = 2:n_sec
    d_edges(k) = sqrt(d_model(k-1) * d_model(k));
end
dw_model = diff(d_edges);

phi_bc      = bc.phi_bc_daily(bc.id_model_best, :);
spec_model  = phi_bc(:) ./ dw_model(:);

d_small = d_model(bc.mask_small);       % model bins < 100 um
s_small = spec_model(bc.mask_small);    % filled by power law

d_vis   = d_model(bc.mask_uvp_model);   % model bins 100-2000 um
s_vis   = spec_model(bc.mask_uvp_model); % mapped from UVP

% ---------------------------------------------------------------
% 4. Figure
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

figure('Units','centimeters','Position',[2 2 12 9]);
hold on;

% power-law line
h1 = plot(d_line, spec_line, 'k--', 'LineWidth', 1);

% UVP data
ok = spec_uvp > 0;
h2 = plot(d_uvp(ok), spec_uvp(ok), 'o', ...
    'Color', [0.2 0.4 0.8], 'MarkerSize', 5, ...
    'MarkerFaceColor', [0.2 0.4 0.8]);

% model bins mapped from UVP (100-2000 um)
h3 = plot(d_vis, s_vis, 's', ...
    'Color', [0.2 0.4 0.8], 'MarkerSize', 5, 'LineWidth', 1.2);

% model bins filled by power law (< 100 um)
ok_s = s_small > 0;
h4 = plot(d_small(ok_s), s_small(ok_s), 's', ...
    'Color', [0.8 0.2 0.2], 'MarkerSize', 5, ...
    'MarkerFaceColor', [0.8 0.2 0.2]);

% mark UVP detection limit
xline(100, ':', 'Color', [0.5 0.5 0.5]);

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([1 3000]);
xlabel('Diameter (\mum)');
ylabel('\phi / \Deltad  (ppmV  \mum^{-1})');
legend([h2 h3 h4 h1], ...
    {sprintf('UVP 100–2000 \\mum (slope %.2f)', slope), ...
     'model bins (UVP mapped)', ...
     'model bins (power-law fill)', ...
     'power-law fit'}, ...
    'location', 'southwest', 'FontSize', 7);
title(sprintf('Volume spectrum at 100 m  —  %d', bc.best_date));

saveas(gcf, fullfile(fig_dir, 'uvp_spectrum_100m.png'));
fprintf('Saved uvp_spectrum_100m.png\n');

% ---------------------------------------------------------------
% 5. Also plot spectra at multiple comparison depths (125-475 m)
%    Shows how the spectrum changes with depth
% ---------------------------------------------------------------
z_cmp  = bc.z_compare;
n_cmp  = numel(z_cmp);
colors = parula(n_cmp);

figure('Units','centimeters','Position',[2 2 12 9]);
hold on;

for i = 1:n_cmp
    [~, iz_u] = min(abs(uvpd.depth_m - z_cmp(i)));
    phi_row = squeeze(uvpd.phi(bc.id_uvp_best, iz_u, mask_uvp));
    phi_row(isnan(phi_row)) = 0;
    spec_row = phi_row(:) ./ dw_uvp(:);
    ok_row = spec_row > 0;
    if any(ok_row)
        plot(d_uvp(ok_row), spec_row(ok_row), 'o-', ...
            'Color', colors(i,:), 'MarkerSize', 3, 'LineWidth', 1, ...
            'DisplayName', sprintf('%d m', round(z_cmp(i))));
    end
end

% reference power-law line (from 100 m fit)
plot(d_line, spec_line, 'k--', 'LineWidth', 1, 'DisplayName', 'fit at 100 m');
xline(100, ':', 'Color', [0.5 0.5 0.5]);

set(gca, 'XScale', 'log', 'YScale', 'log');
xlim([50 3000]);
xlabel('Diameter (\mum)');
ylabel('\phi / \Deltad  (ppmV  \mum^{-1})');
legend('location', 'southwest', 'FontSize', 7);
title(sprintf('UVP spectra by depth  —  %d', bc.best_date));

saveas(gcf, fullfile(fig_dir, 'uvp_spectrum_depth.png'));
fprintf('Saved uvp_spectrum_depth.png\n');

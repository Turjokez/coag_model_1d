% run_sinking_velocity_diagnostic.m
%
% Print sinking velocity w(d) for all 30 model bins using kriest_8 law.
% No simulation needed — just grid + sinking law.
%
% Purpose: check if fast sinking at depth explains low deep standing stock.
% If w > 100 m/day for large bins in the UVP range, particles transit
% through 350-500 m quickly and the steady-state BV will be low.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

cfg = SimulationConfig();
cfg.n_sections  = 30;
cfg.sinking_law = 'kriest_8';
cfg.r_to_rg     = 1.6;
cfg.validate();

grid_c = cfg.derive();

% bin diameters [um]
av_vol = grid_c.av_vol(:);
r_cm   = (0.75 / pi * av_vol).^(1/3);
d_um   = 2 * r_cm * 1e4;
d_cm   = d_um / 1e4;

% sinking velocity [cm/s] -> [m/day]
v_cms  = SettlingVelocityService.velocityForSections(grid_c, cfg);
w_mday = v_cms * cfg.day_to_sec / 100;

% UVP range mask
mask_uvp = d_um >= 100 & d_um < 2000;

fprintf('\nBin  d (um)   w (m/day)   in UVP range\n');
for i = 1:cfg.n_sections
    flag = '';
    if mask_uvp(i), flag = ' <--'; end
    fprintf(' %2d   %6.1f    %7.1f  %s\n', i, d_um(i), w_mday(i), flag);
end

fprintf('\nUVP-range summary (100-2000 um):\n');
fprintf('  min w = %.1f m/day  (bin %d, d=%.0f um)\n', ...
    min(w_mday(mask_uvp)), find(mask_uvp,1,'first'), min(d_um(mask_uvp)));
fprintf('  max w = %.1f m/day  (bin %d, d=%.0f um)\n', ...
    max(w_mday(mask_uvp)), find(mask_uvp,1,'last'), max(d_um(mask_uvp)));

% transit time through 350 m (325->475 m band)
dz_band = 150;   % m
fprintf('\nTransit time through %d m band:\n', dz_band);
fprintf('  smallest UVP bin: %.1f days\n', dz_band / min(w_mday(mask_uvp)));
fprintf('  largest  UVP bin: %.2f days\n', dz_band / max(w_mday(mask_uvp)));

% ---------------------------------------------------------------
% Plot
% ---------------------------------------------------------------
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

figure('Units', 'centimeters', 'Position', [2 2 9 8], 'Color', 'white');
loglog(d_um, w_mday, 'k-o', 'MarkerSize', 3, 'LineWidth', 1.2); hold on;
loglog(d_um(mask_uvp), w_mday(mask_uvp), 'r-o', 'MarkerSize', 3, 'LineWidth', 1.2);
xline(100,  'k--', 'LineWidth', 0.8);
xline(2000, 'k--', 'LineWidth', 0.8);
xlabel('Diameter (\mum)');
ylabel('w (m day^{-1})');
legend('all bins', 'UVP range', 'Location', 'northwest', 'FontSize', 7);
title('Sinking velocity (kriest\_8)', 'FontWeight', 'normal');

saveas(gcf, fullfile(fig_dir, 'sinking_velocity_diagnostic.png'));
fprintf('\nSaved sinking_velocity_diagnostic.png\n');

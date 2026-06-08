% check_uvp_size_at_depth.m
% Check which size bins dominate UVP phi at surface vs depth.
%
% If large bins (>1 mm) dominate at depth -> probably zooplankton,
% not marine snow. Model can't match those.
% If small bins (57-300 um) dominate -> model has a real production gap.

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

uvp_file = fullfile(script_dir, '..', '..', 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

uvp = parse_uvp(uvp_file);

% depth zones
surf_mask  = uvp.depth_m <= 50;
mid_mask   = uvp.depth_m > 50  & uvp.depth_m <= 200;
deep_mask  = uvp.depth_m > 200;

phi_surf  = sum(uvp.phi(surf_mask,  :), 1, 'omitnan');   % 1 x 27
phi_mid   = sum(uvp.phi(mid_mask,   :), 1, 'omitnan');
phi_deep  = sum(uvp.phi(deep_mask,  :), 1, 'omitnan');

% normalize each zone to total so shape is clear
phi_surf_n = phi_surf  / max(sum(phi_surf),  1e-30);
phi_mid_n  = phi_mid   / max(sum(phi_mid),   1e-30);
phi_deep_n = phi_deep  / max(sum(phi_deep),  1e-30);

fprintf('=== UVP size composition by depth zone ===\n');
% fraction in small (<300 um), medium (300-2000 um), large (>2000 um)
small_mask  = uvp.d_um < 300;
medium_mask = uvp.d_um >= 300 & uvp.d_um < 2000;
large_mask  = uvp.d_um >= 2000;

zones = {'surf (0-50m)', 'mid (50-200m)', 'deep (>200m)'};
phi_each = {phi_surf, phi_mid, phi_deep};
for iz = 1:3
    p = phi_each{iz};
    tot = sum(p);
    if tot == 0, continue; end
    fprintf('%s: small=%.1f%%  medium=%.1f%%  large=%.1f%%  total=%.2e\n', ...
        zones{iz}, ...
        100*sum(p(small_mask))/tot, ...
        100*sum(p(medium_mask))/tot, ...
        100*sum(p(large_mask))/tot, ...
        tot);
end

% --- figures ---
fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% normalized size spectrum by depth zone
figure;
semilogx(uvp.d_um, phi_surf_n,  'b-',  'DisplayName', '0-50 m');
hold on;
semilogx(uvp.d_um, phi_mid_n,   'g--', 'DisplayName', '50-200 m');
semilogx(uvp.d_um, phi_deep_n,  'r:',  'DisplayName', '>200 m');
hold off;
xlabel('diameter  [\mum]');
ylabel('\phi / \phi_{total}  (normalized)');
legend('location', 'northwest');
title('UVP size spectrum by depth zone (normalized)');
saveas(gcf, fullfile(fig_dir, 'uvp_size_by_depth.png'));

% absolute phi(z) by size class
phi_small_z  = sum(uvp.phi(:, small_mask),  2, 'omitnan');
phi_medium_z = sum(uvp.phi(:, medium_mask), 2, 'omitnan');
phi_large_z  = sum(uvp.phi(:, large_mask),  2, 'omitnan');

figure;
semilogx(phi_small_z,  uvp.depth_m, 'b-',  'DisplayName', '<300 \mum');
hold on;
semilogx(phi_medium_z, uvp.depth_m, 'g--', 'DisplayName', '300-2000 \mum');
semilogx(phi_large_z,  uvp.depth_m, 'r:',  'DisplayName', '>2000 \mum');
hold off;
set(gca, 'YDir', 'reverse');
xlabel('\phi  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('UVP phi(z) by size class');
saveas(gcf, fullfile(fig_dir, 'uvp_phi_z_by_size.png'));

fprintf('\nFigures saved.\n');

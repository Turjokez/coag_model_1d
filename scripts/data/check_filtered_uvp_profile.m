% check_filtered_uvp_profile.m
% Does UVP phi increase with depth after <2000 um filter?
% Print numbers and plot to check if the trend is real.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));
addpath(fullfile(repo_root, 'scripts', 'data'));

uvp_file = fullfile(repo_root, 'data', 'NA', 'uvp', 'raw', ...
    '745f00117e_EXPORTS-EXPORTSNA_UVP5-ParticulateLevel2_differential_survey_20210504-20210529_R1.sb');

uvp = parse_uvp(uvp_file);

d_max_um = 2000;
mask_agg = uvp.d_um < d_max_um;

% total phi: no filter
phi_all   = uvp.phi;
phi_all(isnan(phi_all)) = 0;
phi_total_all = sum(phi_all, 2);

% aggregate phi: <2000 um only
phi_filt  = uvp.phi(:, mask_agg);
phi_filt(isnan(phi_filt)) = 0;
phi_total_filt = sum(phi_filt, 2);

% fraction that is aggregate at each depth
frac_agg = phi_total_filt ./ max(phi_total_all, 1e-30);

fprintf('depth(m)   phi_all      phi_<2mm     frac_agg\n');
fprintf('--------------------------------------------------\n');
for i = 1:numel(uvp.depth_m)
    fprintf('  %6.1f   %.3e   %.3e   %.2f\n', ...
        uvp.depth_m(i), phi_total_all(i), phi_total_filt(i), frac_agg(i));
end

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end

% two panels: left = phi(z), right = aggregate fraction
figure;

subplot(1,2,1);
semilogy(phi_total_all,  uvp.depth_m, 'r--', 'DisplayName', 'all sizes');
hold on;
semilogy(phi_total_filt, uvp.depth_m, 'b-',  'DisplayName', '<2000 \mum');
hold off;
set(gca, 'YDir', 'reverse');
xlabel('\phi  [cm^3 cm^{-3}]');
ylabel('depth  [m]');
legend('location', 'southeast');
title('UVP phi(z)');

subplot(1,2,2);
plot(frac_agg * 100, uvp.depth_m, 'k-o', 'MarkerSize', 4);
set(gca, 'YDir', 'reverse');
xlabel('aggregate fraction  [%]');
ylabel('depth  [m]');
title('fraction <2000 \mum');

saveas(gcf, fullfile(fig_dir, 'uvp_filtered_profile.png'));
fprintf('\nFigure saved to docs/figures/uvp_filtered_profile.png\n');

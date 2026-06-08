% run_mining_test.m
% Test micro-zoo mining term (Stemmann 2004 Part I, Eq. 25).
%
% Steps:
%   1. Run with mining OFF.
%   2. Run with mining ON.
%   3. Check: mass shifted toward smaller bins, no negatives, fecal produced.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

%% config
cfg = SimulationConfig( ...
    'n_sections',       30, ...
    'sinking_law',      'kriest_8', ...
    'ds_kernel_mode',   'sinking_law', ...
    'enable_coag',      false, ...
    'enable_sinking',   true, ...
    'enable_disagg',    false, ...
    'enable_zoo',       true, ...
    'zoo_Zc',           0, ...
    'zoo_Zf',           0, ...
    'zoo_p',            0.3, ...
    'zoo_ic',           7, ...
    'fp_alpha_cross',   0, ...
    'enable_surface_pp', false, ...
    't_final',          10, ...
    'delta_t',          0.4);

n_z   = 20;
n_sec = cfg.n_sections;

grid  = cfg.derive();
cgrid = ColumnGrid(1000, n_z);
base_prof = DepthProfile.typical(cgrid.z_centers);
prof = DepthProfile(base_prof.z, base_prof.T_K, base_prof.S, ...
    base_prof.rho, base_prof.nu, base_prof.eps, base_prof.Kz);
prof.Zm = 250 * ones(n_z, 1);

%% initial condition: one aggregate size, no production
Y0      = zeros(n_z, n_sec);
Y0(1,15) = 1e-3;

%% run 1: mining OFF
cfg.enable_mining = false;
rhs_off = ColumnRHS(cfg, grid, cgrid, prof);

Y    = Y0;
Yfp  = zeros(n_z, n_sec);
t    = 0;
dt   = cfg.delta_t;
t_end = cfg.t_final;
while t < t_end
    [Y, Yfp] = rhs_off.stepY(Y, dt, Yfp);
    t = t + dt;
end
Y_off  = Y;
Yfp_off = Yfp;

%% run 2: mining ON
cfg.enable_mining = true;
rhs_on = ColumnRHS(cfg, grid, cgrid, prof);

Y    = Y0;
Yfp  = zeros(n_z, n_sec);
t    = 0;
while t < t_end
    [Y, Yfp] = rhs_on.stepY(Y, dt, Yfp);
    t = t + dt;
end
Y_on   = Y;
Yfp_on = Yfp;

%% checks
assert(all(Y_on(:) >= 0),  'FAIL: negatives in Y with mining on');
assert(all(Yfp_on(:) >= 0),'FAIL: negatives in Yfp with mining on');

total_off = sum(Y_off(:)) + sum(Yfp_off(:));
total_on  = sum(Y_on(:))  + sum(Yfp_on(:));
fprintf('Total mass (mining off): %.4e\n', total_off);
fprintf('Total mass (mining on):  %.4e\n', total_on);
fprintf('Mining reduces total by: %.1f%%\n', 100*(total_off - total_on)/total_off);

% size spectrum: mining should shift mass toward smaller bins
spec_off = sum(Y_off, 1);   % total across depths
spec_on  = sum(Y_on,  1);
[~, peak_off] = max(spec_off);
[~, peak_on]  = max(spec_on);
fprintf('Peak bin (mining off): %d\n', peak_off);
fprintf('Peak bin (mining on):  %d\n', peak_on);
fprintf('Fecal produced:        %.4e\n', sum(Yfp_on(:)));

assert(total_on < total_off, 'FAIL: mining did not reduce total mass');
assert(total_on > 0.5 * total_off, 'FAIL: mining is too strong in this test');
assert(peak_on <= peak_off, 'FAIL: mining did not shift spectrum smaller');
assert(sum(Yfp_on(:)) > sum(Yfp_off(:)), 'FAIL: mining did not produce fecal material');

%% figure: size spectrum at t=60 days
figure;
bins = 1:n_sec;
plot(bins, spec_off, 'k-',  'DisplayName', 'mining off'); hold on;
plot(bins, spec_on,  'b--', 'DisplayName', 'mining on');
xlabel('bin');
ylabel('total biovolume');
legend;
title('size spectrum at t=10 d');

fig_dir = fullfile(repo_root, 'docs', 'figures');
if ~exist(fig_dir, 'dir'), mkdir(fig_dir); end
saveas(gcf, fullfile(fig_dir, 'mining_test_spectrum.png'));

fprintf('\nAll checks passed.\n');

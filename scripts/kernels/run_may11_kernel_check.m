% run_may11_kernel_check
% Check a few kernel values by hand and compare with model matrices.

clear; close all; clc;

addpath('src');
repo_root = pwd;
if ~exist('SimulationConfig', 'class')
    repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    addpath(genpath(fullfile(repo_root, 'src')));
end

% Build config and size grid.
cfg = SimulationConfig('n_sections', 20, 'sinking_law', 'kriest_8', ...
    'ds_kernel_mode', 'sinking_law');
grid = cfg.derive();

% Section pairs to check.
pairs = [1 5;
         1 10;
         5 15];

% Diameters [cm].
if isprop(grid, 'd_low') && ~isempty(grid.d_low)
    d_sec = grid.d_low(:);
else
    d_sec = grid.getVolumeDiameters();
end

% Physical constants (same as config defaults).
kT    = 1.3e-16 * (20 + 273);   % erg
mu    = 1.0275 * 0.01;          % g cm^-1 s^-1
gamma = 0.1;                    % s^-1

% Settling speeds [cm/s], one per section.
w_sec = SettlingVelocityService.velocityForSections(grid, cfg);

fprintf('=== Hand calculation (cm^3/s) ===\n');
for p = 1:size(pairs, 1)
    i = pairs(p, 1);
    j = pairs(p, 2);
    d_i = d_sec(i);
    d_j = d_sec(j);

    beta_brown = (2 * kT) / (3 * mu) * ((d_i + d_j)^2 / (d_i * d_j));
    beta_shear = (pi / 6) * gamma * ((d_i + d_j) / 2)^3;
    beta_ds    = (pi / 4) * (d_i + d_j)^2 * abs(w_sec(i) - w_sec(j));

    fprintf('pair (%d,%d):  d_i=%.2e cm  d_j=%.2e cm\n', i, j, d_i, d_j);
    fprintf('  beta_brown = %.3e cm3/s\n', beta_brown);
    fprintf('  beta_shear = %.3e cm3/s\n', beta_shear);
    fprintf('  beta_ds    = %.3e cm3/s\n', beta_ds);
    fprintf('  DS/Brown   = %.2f\n', beta_ds / beta_brown);
    fprintf('  DS/Shear   = %.2f\n', beta_ds / beta_shear);
    fprintf('\n');
end

% Model kernels from assembler.
assembler = BetaAssembler(cfg, grid);
b_brown = assembler.computeFor('KernelBrown');
b_shear = assembler.computeFor('KernelCurSh');
b_ds    = assembler.computeFor('KernelCurDSSinkingLaw');

alpha = 1.0;
day_to_sec = 8.64e4;

fprintf('=== Model b25 values and scaled comparison ===\n');
for p = 1:size(pairs, 1)
    i = pairs(p, 1);
    j = pairs(p, 2);
    d_i = d_sec(i);
    d_j = d_sec(j);

    beta_brown = (2 * kT) / (3 * mu) * ((d_i + d_j)^2 / (d_i * d_j));
    beta_shear = (pi / 6) * gamma * ((d_i + d_j) / 2)^3;
    beta_ds    = (pi / 4) * (d_i + d_j)^2 * abs(w_sec(i) - w_sec(j));

    % raw b25 from model
    mb_raw = b_brown.b25(i, j);
    ms_raw = b_shear.b25(i, j);
    md_raw = b_ds.b25(i, j);

    % add the same scaling used by model assembly for each mechanism
    mb_scaled = mb_raw * grid.conBr * alpha * day_to_sec;
    ms_scaled = ms_raw * cfg.gamma * alpha * day_to_sec;
    md_scaled = md_raw * grid.setcon * alpha * day_to_sec;

    rb = mb_scaled / (beta_brown * alpha * day_to_sec);
    rs = ms_scaled / (beta_shear * alpha * day_to_sec);
    rd = md_scaled / (beta_ds * alpha * day_to_sec);

    fprintf('pair (%d,%d)\n', i, j);
    fprintf('  model b25 brown raw = %.3e\n', mb_raw);
    fprintf('  model b25 shear raw = %.3e\n', ms_raw);
    fprintf('  model b25 ds raw    = %.3e\n', md_raw);
    fprintf('  ratio brown model/(hand*alpha*day) = %.3f\n', rb);
    fprintf('  ratio shear model/(hand*alpha*day) = %.3f\n', rs);
    fprintf('  ratio ds model/(hand*alpha*day)    = %.3f\n', rd);
    fprintf('\n');
end


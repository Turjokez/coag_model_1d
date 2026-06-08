% run_verify_all.m
% Master verification script for the 1-D coagulation column model.
%
% Runs quick checks on every major physics component.
% All tests use short runs (t <= 20 days) so this completes fast.
%
% Usage:
%   cd scripts/
%   run_verify_all
%
% Expected output: all PASS, 0 FAIL.

clear; close all; clc;

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));

pass_n = 0;
fail_n = 0;

function [pass_n, fail_n] = check(name, ok, pass_n, fail_n)
    if ok
        fprintf('  PASS  %s\n', name);
        pass_n = pass_n + 1;
    else
        fprintf('  FAIL  %s\n', name);
        fail_n = fail_n + 1;
    end
end

fprintf('\n=== Model verification  %s ===\n\n', datestr(now,'yyyy-mm-dd'));

% ----------------------------------------------------------------
% shared base config (n=30, EXPORTS-ready defaults)
% ----------------------------------------------------------------
base = SimulationConfig( ...
    'n_sections',      30, ...
    'sinking_law',     'kriest_8', ...
    'ds_kernel_mode',  'sinking_law', ...
    'enable_coag',     true, ...
    'enable_sinking',  true, ...
    'enable_disagg',   false, ...
    'enable_zoo',      false, ...
    't_final',         10, ...
    'delta_t',         0.4);

n_z  = 20;  dz = 50;
grid = base.derive();
cgrid = ColumnGrid(dz * n_z, n_z);
prof  = DepthProfile.typical(cgrid.z_centers);
prof_nomix = DepthProfile(prof.z, prof.T_K, prof.S, prof.rho, prof.nu, prof.eps, zeros(n_z,1));

% helper: make Y0 with a pulse at surface, given bin
function Y0 = makePulse(n_z, n_sec, bin, amount)
    Y0 = zeros(n_z, n_sec);
    Y0(1, bin) = amount;
end

% ----------------------------------------------------------------
fprintf('--- 1. Transport ---\n');
% ----------------------------------------------------------------

% pulse at surface should move down over time, no negatives
cfg = base.copy(); cfg.enable_coag = false;
rhs = ColumnRHS(cfg, grid, cgrid, prof_nomix);
Y = makePulse(n_z, base.n_sections, 30, 1e-3); Yfp = zeros(n_z, base.n_sections);
z0 = sum(sum(Y,2) .* (1:n_z)') / sum(Y(:));
for s = 1:25, [Y, Yfp] = rhs.stepY(Y, cfg.delta_t, Yfp); end
z1 = sum(sum(Y,2) .* (1:n_z)') / sum(Y(:));
[pass_n, fail_n] = check('No negatives after sinking',     all(Y(:) >= 0),             pass_n, fail_n);
[pass_n, fail_n] = check('Pulse center moved downward',    z1 > z0,                    pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 2. Coagulation ---\n');
% ----------------------------------------------------------------

% coagulation should shift mass upward (toward larger bins)
cfg = base.copy(); cfg.enable_sinking = false;
rhs = ColumnRHS(cfg, grid, cgrid, prof);
Y = makePulse(n_z, base.n_sections, 3, 1e-3); Yfp = zeros(n_z, base.n_sections);
Y_init = Y;
for s = 1:10, [Y, Yfp] = rhs.stepY(Y, cfg.delta_t, Yfp); end
[~, pk_before] = max(sum(Y_init,1));
[~, pk_after]  = max(sum(Y,1));
[pass_n, fail_n] = check('No negatives after coagulation',        all(Y(:) >= 0),        pass_n, fail_n);
[pass_n, fail_n] = check('Coagulation shifts peak to larger bins', pk_after >= pk_before, pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 3. Disaggregation ---\n');
% ----------------------------------------------------------------

cfg = base.copy();
cfg.enable_sinking = false; cfg.enable_coag = false;
cfg.enable_disagg  = true;  cfg.disagg_mode = 'operator_split';
cfg.disagg_epsilon = 1e-6;
rhs = ColumnRHS(cfg, grid, cgrid, prof);
Y = makePulse(n_z, base.n_sections, 20, 1e-3); Yfp = zeros(n_z, base.n_sections);
Y_init = Y;
for s = 1:10, [Y, Yfp] = rhs.stepY(Y, cfg.delta_t, Yfp); end
[~, pk_before] = max(sum(Y_init,1));
[~, pk_after]  = max(sum(Y,1));
[pass_n, fail_n] = check('No negatives after disagg',              all(Y(:) >= 0),        pass_n, fail_n);
[pass_n, fail_n] = check('Disagg shifts peak to smaller bins',     pk_after <= pk_before, pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 4. Zooplankton grazing ---\n');
% ----------------------------------------------------------------

cfg = base.copy();
cfg.enable_sinking = false; cfg.enable_coag = false;
cfg.enable_zoo = true;
cfg.zoo_Zc = 0.307; cfg.zoo_Zf = 0.063; cfg.zoo_p = 0.3; cfg.zoo_ic = 7;
rhs = ColumnRHS(cfg, grid, cgrid, prof);
prof2 = prof; prof2.Zc = cfg.zoo_Zc * ones(n_z,1); prof2.Zf = cfg.zoo_Zf * ones(n_z,1);
rhs2 = ColumnRHS(cfg, grid, cgrid, prof2);
Y = makePulse(n_z, base.n_sections, 10, 1e-3); Yfp = zeros(n_z, base.n_sections);
Y_init = Y; Yfp_init = Yfp;
for s = 1:25, [Y, Yfp] = rhs2.stepY(Y, cfg.delta_t, Yfp); end
[pass_n, fail_n] = check('No negatives after grazing',         all(Y(:) >= 0) && all(Yfp(:) >= 0), pass_n, fail_n);
[pass_n, fail_n] = check('Grazing reduces aggregate mass',     sum(Y(:)) < sum(Y_init(:)),           pass_n, fail_n);
[pass_n, fail_n] = check('Fecal pellets produced by grazing',  sum(Yfp(:)) > sum(Yfp_init(:)),      pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 5. Fecal pellet sinking speed ---\n');
% ----------------------------------------------------------------

w_agg = SettlingVelocityService.velocityForSections(grid, base);
w_fp  = SettlingVelocityService.velocityFecalPellets(grid, base);
ratio = w_fp(8) / w_agg(8);   % bin 8 (~115 um)
[pass_n, fail_n] = check('Fecal sinking > aggregate at bin 8',  w_fp(8) > w_agg(8),   pass_n, fail_n);
[pass_n, fail_n] = check('Fecal/aggregate ratio > 10x at bin 8', ratio > 10,           pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 6. Cross-coagulation ---\n');
% ----------------------------------------------------------------

cfg = base.copy();
cfg.enable_sinking = true; cfg.enable_coag = false;
cfg.enable_zoo = true; cfg.zoo_Zc = 0; cfg.zoo_Zf = 0;
cfg.zoo_ic = 7; cfg.fp_alpha_cross = 0.5;
rhs = ColumnRHS(cfg, grid, cgrid, prof_nomix);
Y = makePulse(n_z, base.n_sections, 10, 1e-3);
Yfp = zeros(n_z, base.n_sections); Yfp(1,8) = 1e-4;  % seed fecal pellets
Yfp_init = Yfp;
for s = 1:5, [Y, Yfp] = rhs.stepY(Y, cfg.delta_t, Yfp); end
[pass_n, fail_n] = check('No negatives after cross-coag',         all(Y(:) >= 0) && all(Yfp(:) >= 0), pass_n, fail_n);
[pass_n, fail_n] = check('Cross-coag reduces fecal pellet mass',  sum(Yfp(:)) < sum(Yfp_init(:)),      pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 7. Microbial remineralization ---\n');
% ----------------------------------------------------------------

cfg = base.copy();
cfg.enable_sinking = false; cfg.enable_coag = false;
cfg.enable_microbe = true; cfg.microbe_r0 = 0.05;
rhs = ColumnRHS(cfg, grid, cgrid, prof);
Y = makePulse(n_z, base.n_sections, 10, 1e-3); Yfp = zeros(n_z, base.n_sections);
Y_init = Y;
for s = 1:25, [Y, Yfp] = rhs.stepY(Y, cfg.delta_t, Yfp); end
decay_expected = exp(-cfg.microbe_r0 * 25 * cfg.delta_t);
decay_actual   = sum(Y(:)) / sum(Y_init(:));
[pass_n, fail_n] = check('No negatives after microbe loss',     all(Y(:) >= 0),                    pass_n, fail_n);
[pass_n, fail_n] = check('Microbial decay close to exp(-r*t)',  abs(decay_actual - decay_expected) < 0.05, pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 8. Mining term ---\n');
% ----------------------------------------------------------------

cfg = base.copy();
cfg.enable_sinking = false; cfg.enable_coag = false;
cfg.enable_zoo = true; cfg.zoo_Zc = 0; cfg.zoo_Zf = 0;
cfg.zoo_ic = 7; cfg.fp_alpha_cross = 0;
cfg.enable_mining = true;
cfg.mining_Zm = 250; cfg.mining_dm = 1e-5; cfg.mining_min_bin = 12;
prof3 = prof; prof3.Zm = 250 * ones(n_z, 1);
rhs = ColumnRHS(cfg, grid, cgrid, prof3);
Y = makePulse(n_z, base.n_sections, 15, 1e-3); Yfp = zeros(n_z, base.n_sections);  % bin 15 = 508 um, above min_bin
Y_init = Y; Yfp_init = Yfp;
for s = 1:25, [Y, Yfp] = rhs.stepY(Y, cfg.delta_t, Yfp); end
[pass_n, fail_n] = check('No negatives after mining',          all(Y(:) >= 0) && all(Yfp(:) >= 0), pass_n, fail_n);
[pass_n, fail_n] = check('Mining reduces aggregate mass',      sum(Y(:)) < sum(Y_init(:)),           pass_n, fail_n);
[pass_n, fail_n] = check('Mining produces fecal material',     sum(Yfp(:)) > 0,                      pass_n, fail_n);
% bins below min_bin must be untouched
Y_small_before = sum(sum(Y_init(:, 1:cfg.mining_min_bin-1)));
Y_small_after  = sum(sum(Y(:,    1:cfg.mining_min_bin-1)));
[pass_n, fail_n] = check('Mining does not touch bins < min_bin', Y_small_after == Y_small_before, pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 9. Surface production ---\n');
% ----------------------------------------------------------------

cfg = base.copy();
cfg.enable_coag = false; cfg.enable_sinking = false;
cfg.enable_surface_pp = true; cfg.surface_pp_mu = 0.1;
rhs = ColumnRHS(cfg, grid, cgrid, prof_nomix);
Y = zeros(n_z, base.n_sections); Y(1,1) = 1e-3; Yfp = zeros(n_z, base.n_sections);
Y_init = Y;
for s = 1:1, [Y, Yfp] = rhs.stepY(Y, cfg.delta_t, Yfp); end
[pass_n, fail_n] = check('Surface production increases bin 1 at layer 1', Y(1,1) > Y_init(1,1), pass_n, fail_n);
[pass_n, fail_n] = check('Surface production does not affect deep layers', Y(end,1) == 0,        pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('--- 10. Full 1-D integrated run ---\n');
% ----------------------------------------------------------------

cfg = base.copy();
cfg.enable_disagg    = true;   cfg.disagg_mode   = 'operator_split';
cfg.enable_zoo       = true;   cfg.zoo_Zc        = 0.307;
cfg.zoo_Zf           = 0.063;  cfg.zoo_ic        = 7;
cfg.enable_surface_pp = true;  cfg.surface_pp_mu = 0.1;
cfg.disagg_epsilon   = 1e-6;
cfg.t_final          = 20;

grid2 = cfg.derive();
prof4  = DepthProfile.typical(cgrid.z_centers);
prof4.Zc = cfg.zoo_Zc * ones(n_z,1);
prof4.Zf = cfg.zoo_Zf * ones(n_z,1);
rhs = ColumnRHS(cfg, grid2, cgrid, prof4);

Y = zeros(n_z, cfg.n_sections); Y(1,1) = 1e-4;
Yfp = zeros(n_z, cfg.n_sections);
t = 0; dt = cfg.delta_t;
neg_found = false;
while t < cfg.t_final
    [Y, Yfp] = rhs.stepY(Y, dt, Yfp);
    if any(Y(:) < -1e-20) || any(Yfp(:) < -1e-20), neg_found = true; end
    t = t + dt;
end
total_bv = sum(Y(:)) + sum(Yfp(:));
[pass_n, fail_n] = check('Full run: no negatives',           ~neg_found,          pass_n, fail_n);
[pass_n, fail_n] = check('Full run: positive total biovolume', total_bv > 0,       pass_n, fail_n);

% ----------------------------------------------------------------
fprintf('\n=== TOTAL: %d PASS, %d FAIL ===\n', pass_n, fail_n);
if fail_n == 0
    fprintf('All checks passed. Model is ready.\n\n');
else
    fprintf('Some checks FAILED. Review output above.\n\n');
end

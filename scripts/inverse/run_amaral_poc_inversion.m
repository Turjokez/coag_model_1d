% run_amaral_poc_inversion.m
%
% Two-class POC cycling inversion (Amaral et al. 2022) for EXPORTS-NA.
%
% Model (steady state, per layer):
%   wS*(PS_bot - PS_top) = [-rem_S*PS_bar - agg*PS_bar^2 + dis*PL_bar]*dz + eS*dz
%   wL*(PL_bot - PL_top) = [ agg*PS_bar^2 - rem_L*PL_bar - dis*PL_bar]*dz + eL*dz
%
% eS, eL absorb unresolved processes (DVM, lateral transport, etc.).
%
% ATI cost:
%   J = sum[(p - p0)^2 / sig0^2] + sum[(eS_k/sig_e)^2 + (eL_k/sig_e)^2]
% where eS, eL are computed from the equations given p (not free variables).
%
% Units: POC mmol/m^3, Th-234 dpm/L (-> ratio x1000 = dpm/mmol),
%        w m/day, rates day^-1, agg m^3/(mmol*day)

clear; close all;

pump_file = fullfile('..', '..', 'data', 'NA', 'thorium_buesseler', 'raw', ...
    'a026c22a81_EXPORTS_EXPORTSNA_DY131_IN_SITU_PUMPS_SURVEY_R1.sb');

% ── 1. Parse pump file ──────────────────────────────────────────────────────
[fields, units, data] = parse_sb(pump_file);

col = @(name) find(strcmp(fields, name));

depth_all = data(:, col('depth'));

% POC (mmol/m^3)
PC1 = data(:, col('PC_1umfilt_5umprefilt'));
PC2 = data(:, col('PC_5umfilt_51umprefilt'));
PC3 = data(:, col('PC_51umfilt_335umprefilt'));
PC4 = data(:, col('PC_335umfilt'));
PS_all = PC1 + PC2;   % small, 1-51 um
PL_all = PC3 + PC4;   % large, >51 um

% Th-234 (dpm/L)
T1 = data(:, col('conc_Th_234_1umfilt_5umprefilt'));
T2 = data(:, col('conc_Th_234_5umfilt_51umprefilt'));
T3 = data(:, col('conc_Th_234_51umfilt_335umprefilt'));
T4 = data(:, col('conc_Th_234_335umfilt'));
Th_all = T1 + T2 + T3 + T4;   % total particulate, dpm/L

% ── 2. Cruise-mean profiles ─────────────────────────────────────────────────
z_target = [20 50 75 95 125 175 330 500]';
PS_mean  = depth_mean(depth_all, PS_all, z_target);
PL_mean  = depth_mean(depth_all, PL_all, z_target);
Th_mean  = depth_mean(depth_all, Th_all, z_target);

fprintf('Cruise-mean profiles:\n');
fprintf('%6s  %7s  %7s  %10s\n', 'z(m)', 'PS', 'PL', 'Th(dpm/L)');
for i = 1:length(z_target)
    fprintf('  %4.0f  %7.3f  %7.3f  %10.4f\n', ...
        z_target(i), PS_mean(i), PL_mean(i), Th_mean(i));
end

% ── 3. Mesopelagic layers: 95, 125, 175, 330, 500 m ────────────────────────
z_idx = [4 5 6 7 8];   % indices for 95,125,175,330,500
z_bnd = z_target(z_idx);
PS_bnd = PS_mean(z_idx);
PL_bnd = PL_mean(z_idx);
Th_bnd = Th_mean(z_idx);

dz     = diff(z_bnd);            % layer thicknesses [m]
PS_bar = 0.5*(PS_bnd(1:end-1) + PS_bnd(2:end));
PL_bar = 0.5*(PL_bnd(1:end-1) + PL_bnd(2:end));
n_lyr  = length(dz);

fprintf('\n%d mesopelagic layers: %s m\n', n_lyr, ...
    strjoin(arrayfun(@(a,b) sprintf('%d-%dm',a,b), z_bnd(1:end-1)', z_bnd(2:end)', ...
    'UniformOutput', false), ', '));

% ── 4. ATI cost function ────────────────────────────────────────────────────
% Prior (Amaral Table 2, NP station).
% agg: 0.003 dm^3/(mmol*day) = 3e-6 m^3/(mmol*day)
p0   = [2.0;  20.0;  0.10;  0.15;  3e-6;  0.43];
p_sig = [1.0;  10.0;  0.05;  0.075; 3e-7;  0.05];

% sigma_e: allowed unresolved flux divergence (mmol/m^3/day)
sig_e = 0.30;

% Cost function handle
cost_fn = @(p) ati_cost(p, PS_bnd, PL_bnd, PS_bar, PL_bar, dz, p0, p_sig, sig_e);

% ── 5. Optimize ─────────────────────────────────────────────────────────────
% Multi-start from prior and perturbed starts
lb = [0.1; 1.0; 0.0; 0.0; 0.0;  0.0];
ub = [10;  200; 2.0; 2.0; 1e-4; 5.0];

opts = optimoptions('fmincon', 'Display', 'off', ...
    'MaxFunctionEvaluations', 2e4, 'OptimalityTolerance', 1e-12);

p_best = p0;
J_best = inf;

% 10 random starts
rng(42);
starts = [p0, p0 .* (0.5 + rand(6, 9))];
for k = 1:size(starts, 2)
    p_try = max(lb, min(ub, starts(:, k)));
    try
        [p_out, J_out] = fmincon(cost_fn, p_try, [], [], [], [], lb, ub, [], opts);
        if J_out < J_best
            J_best = J_out;
            p_best = p_out;
        end
    catch
    end
end

p_fit = p_best;
wS = p_fit(1); wL = p_fit(2);
rem_S = p_fit(3); rem_L = p_fit(4);
agg = p_fit(5); dis = p_fit(6);
p_names = {'wS', 'wL', 'rem_S', 'rem_L', 'agg', 'dis'};
p_units = {'m/day', 'm/day', 'day^-1', 'day^-1', 'm3/(mmol*day)', 'day^-1'};

fprintf('\n── Fitted parameters vs prior ─────────────────────────────────\n');
for i = 1:6
    pull = (p_fit(i) - p0(i)) / p_sig(i);
    fprintf('  %-7s: %10.4g  prior=%g +/- %g  pull=%+.1fsig  [%s]\n', ...
        p_names{i}, p_fit(i), p0(i), p_sig(i), pull, p_units{i});
end

% Decompose cost
[eS, eL] = model_errors(p_fit, PS_bnd, PL_bnd, PS_bar, PL_bar, dz);
J_param = sum(((p_fit - p0)./p_sig).^2);
J_eq    = sum((eS/sig_e).^2) + sum((eL/sig_e).^2);
fprintf('  Cost J = %.4f  (J_param=%.2f, J_eq=%.2f)\n', J_best, J_param, J_eq);

fprintf('\n── Model errors (eS, eL) per layer ────────────────────────────\n');
fprintf('  (eS<0: net sink on small POC not in model; eL>0: unresolved source)\n');
for k = 1:n_lyr
    fprintf('  Layer %d (%d-%dm): eS=%+.3f, eL=%+.3f  mmol/m3/day\n', ...
        k, z_bnd(k), z_bnd(k+1), eS(k), eL(k));
end

% ── 6. Layer budgets ────────────────────────────────────────────────────────
fprintf('\n── Small POC budget (mmol/m2/day) ─────────────────────────────\n');
fprintf('  %-14s  %9s  %9s  %8s  %8s  %10s\n', ...
    'Layer', 'Flux_top', 'Flux_bot', 'Rem', 'Dis_in', 'Unresolved');
for k = 1:n_lyr
    fprintf('  %d-%dm:   %9.3f  %9.3f  %8.3f  %8.3f  %10.3f\n', ...
        z_bnd(k), z_bnd(k+1), ...
        wS*PS_bnd(k), wS*PS_bnd(k+1), ...
        rem_S*PS_bar(k)*dz(k), dis*PL_bar(k)*dz(k), eS(k)*dz(k));
end

fprintf('\n── Large POC budget (mmol/m2/day) ─────────────────────────────\n');
fprintf('  %-14s  %9s  %9s  %8s  %8s  %10s\n', ...
    'Layer', 'Flux_top', 'Flux_bot', 'Rem', 'Dis_out', 'Unresolved');
for k = 1:n_lyr
    fprintf('  %d-%dm:   %9.3f  %9.3f  %8.3f  %8.3f  %10.3f\n', ...
        z_bnd(k), z_bnd(k+1), ...
        wL*PL_bnd(k), wL*PL_bnd(k+1), ...
        rem_L*PL_bar(k)*dz(k), dis*PL_bar(k)*dz(k), eL(k)*dz(k));
end

% ── 7. Export flux and Th-234 comparison ───────────────────────────────────
% Th units: dpm/L.  POC units: mmol/m^3.
% Th/POC ratio [dpm/mmol] = Th[dpm/L] / POC[mmol/m^3] * 1000
POC_100 = PS_bnd(1) + PL_bnd(1);           % total POC at 95m
ThPOC   = Th_bnd(1) / POC_100 * 1000;      % dpm/mmol

F_S = wS * PS_bnd(1);
F_L = wL * PL_bnd(1);
F_POC = F_S + F_L;
F_Th_model  = F_POC * ThPOC;
F_Th_bottle = 1098.0;   % dpm/m^2/day from lambda*integral(U-Th) 0-100m
F_POC_Th    = F_Th_bottle / ThPOC;          % Th-implied POC flux

F_POC_500 = wS*PS_bnd(end) + wL*PL_bnd(end);
TE = F_POC_500 / F_POC * 100;

fprintf('\n── Export flux at ~100 m ───────────────────────────────────────\n');
fprintf('  wS*PS(95m) = %.2f * %.3f = %.2f mmol/m2/day\n', wS, PS_bnd(1), F_S);
fprintf('  wL*PL(95m) = %.2f * %.3f = %.2f mmol/m2/day\n', wL, PL_bnd(1), F_L);
fprintf('  Total POC flux (model)   = %.2f mmol/m2/day\n', F_POC);
fprintf('  Th/POC at 95m            = %.1f dpm/mmol\n', ThPOC);
fprintf('  Th flux (model x ratio)  = %.0f dpm/m2/day\n', F_Th_model);
fprintf('  Th flux (bottle deficit) = %.0f dpm/m2/day\n', F_Th_bottle);
fprintf('  ** Ratio model/bottle    = %.1f x **\n', F_Th_model/F_Th_bottle);
fprintf('\n  Th-implied POC flux      = %.2f mmol/m2/day\n', F_POC_Th);
fprintf('  wL needed to match Th    = %.1f m/day\n', ...
    (F_Th_bottle/ThPOC - wS*PS_bnd(1)) / PL_bnd(1));
fprintf('  Transfer efficiency      = %.1f%%  (500/100m)\n', TE);

% ── 8. Size-fraction breakdown and swimmer sensitivity at 95 m ─────────────
% cruise-mean size fractions (from depth_mean on raw columns)
PC_1_5    = depth_mean(depth_all, data(:, col('PC_1umfilt_5umprefilt')),    z_target);
PC_5_51   = depth_mean(depth_all, data(:, col('PC_5umfilt_51umprefilt')),   z_target);
PC_51_335 = depth_mean(depth_all, data(:, col('PC_51umfilt_335umprefilt')), z_target);
PC_335    = depth_mean(depth_all, data(:, col('PC_335umfilt')),             z_target);

idx95 = find(z_target == 95);
PL_51_335 = PC_51_335(idx95);   % 51-335 um (small zooplankton range)
PL_gt335  = PC_335(idx95);      % >335 um (large swimmers)
PS_95_val = PS_bnd(1);

fprintf('\n── Size fractions at 95 m ──────────────────────────────────────\n');
fprintf('  1-5 um:      %.3f mmol/m3\n', PC_1_5(idx95));
fprintf('  5-51 um:     %.3f mmol/m3\n', PC_5_51(idx95));
fprintf('  51-335 um:   %.3f mmol/m3  <- small-zoo range\n', PL_51_335);
fprintf('  >335 um:     %.3f mmol/m3  <- obvious swimmers\n', PL_gt335);
fprintf('  PS (1-51):   %.3f mmol/m3\n', PS_95_val);
fprintf('  PL total:    %.3f mmol/m3\n', PL_bnd(1));

% Flux sensitivity: vary what fraction of PL is treated as sinking
fprintf('\n── Swimmer sensitivity (wS=%.2f, wL=%.2f, Th/POC=%.0f dpm/mmol) ──\n', ...
    wS, wL, ThPOC);
fprintf('  %-35s  %7s  %8s  %10s  %6s\n', ...
    'Scenario', 'PL used', 'F_POC', 'F_Th_model', 'ratio');

scenarios = {
    'Full PL (51-335 + >335)',    PL_bnd(1);
    'Only 51-335 um (remove >335)', PL_51_335;
    'PL x 50%',                   PL_bnd(1)*0.50;
    'PL x 25%',                   PL_bnd(1)*0.25;
    'PL x 10%',                   PL_bnd(1)*0.10;
    'Th-matched PL',              (F_Th_bottle/ThPOC - wS*PS_95_val)/wL;
};
for k = 1:size(scenarios,1)
    PL_use = scenarios{k,2};
    fp = wS*PS_95_val + wL*PL_use;
    ft = fp * ThPOC;
    fprintf('  %-35s  %7.3f  %8.2f  %10.0f  %6.2fx\n', ...
        scenarios{k,1}, PL_use, fp, ft, ft/F_Th_bottle);
end

PL_match = (F_Th_bottle/ThPOC - wS*PS_95_val) / wL;
frac_swimmer = (PL_bnd(1) - PL_match) / PL_bnd(1) * 100;
fprintf('\n  To match Th: PL_sinking = %.3f mmol/m3 (%.0f%% of measured PL)\n', ...
    PL_match, 100 - frac_swimmer);
fprintf('  Implies ~%.0f%% of PL at 95m is non-sinking (swimmers)\n', frac_swimmer);

% ── 9. Simple figure ────────────────────────────────────────────────────────
set(0, 'DefaultAxesFontName', 'Arial', 'DefaultTextFontName', 'Arial');

figure;

% Panel 1: POC profiles
subplot(1, 2, 1);
plot(PS_bnd, z_bnd, 'o-b', 'LineWidth', 1.2, 'MarkerFaceColor', 'b'); hold on;
plot(PL_bnd, z_bnd, 's-r', 'LineWidth', 1.2, 'MarkerFaceColor', 'r');
set(gca, 'YDir', 'reverse');
xlabel('POC (mmol m^{-3})');
ylabel('Depth (m)');
legend('PS (1-51 \mum)', 'PL (>51 \mum)', 'Location', 'southeast');
title('Pump POC');

% Panel 2: model errors
subplot(1, 2, 2);
z_mid = 0.5*(z_bnd(1:end-1) + z_bnd(2:end));
barh(z_mid, eS, 'FaceColor', [0.2 0.5 0.8], 'BarWidth', 0.3); hold on;
barh(z_mid + 40, eL, 'FaceColor', [0.8 0.3 0.3], 'BarWidth', 0.3);
xline(0, 'k--');
set(gca, 'YDir', 'reverse');
xlabel('e (mmol m^{-3} day^{-1})');
ylabel('Depth (m)');
legend('e_S', 'e_L', 'Location', 'southeast');
title('Unresolved terms');

% ── Local functions ─────────────────────────────────────────────────────────
function J = ati_cost(p, PS_bnd, PL_bnd, PS_bar, PL_bar, dz, p0, p_sig, sig_e)
    if any(p < 0), J = 1e9; return; end
    J_param = sum(((p - p0) ./ p_sig).^2);
    [eS, eL] = model_errors(p, PS_bnd, PL_bnd, PS_bar, PL_bar, dz);
    J_eq = sum((eS/sig_e).^2) + sum((eL/sig_e).^2);
    J = J_param + J_eq;
end

function [eS, eL] = model_errors(p, PS_bnd, PL_bnd, PS_bar, PL_bar, dz)
    wS = p(1); wL = p(2);
    rem_S = p(3); rem_L = p(4); agg = p(5); dis = p(6);
    src_S = (-rem_S*PS_bar - agg*PS_bar.^2 + dis*PL_bar) .* dz;
    src_L = ( agg*PS_bar.^2 - rem_L*PL_bar - dis*PL_bar) .* dz;
    eS = (wS*diff(PS_bnd) - src_S) ./ dz;
    eL = (wL*diff(PL_bnd) - src_L) ./ dz;
end

function m = depth_mean(z_all, vals, z_target, tol)
    if nargin < 4, tol = 10; end
    m = nan(size(z_target));
    for i = 1:length(z_target)
        mask = abs(z_all - z_target(i)) <= tol & isfinite(vals);
        if any(mask), m(i) = mean(vals(mask)); end
    end
end

function [fields, units, data] = parse_sb(path)
    fields = {}; units = {}; rows = {};
    fid = fopen(path, 'r');
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        if startsWith(line, '/fields=')
            fields = strsplit(line(9:end), ',');
        elseif startsWith(line, '/units=')
            units = strsplit(line(8:end), ',');
        elseif ~isempty(line) && line(1) ~= '/' && line(1) ~= '!'
            rows{end+1} = strsplit(line, ',');  %#ok<AGROW>
        end
    end
    fclose(fid);
    n_col = length(fields);
    n_row = length(rows);
    data  = nan(n_row, n_col);
    MISS  = [-9999 -8888];
    for r = 1:n_row
        row = rows{r};
        for c = 1:min(n_col, length(row))
            v = str2double(row{c});
            if ~isnan(v) && ~any(v == MISS)
                data(r, c) = v;
            end
        end
    end
end

% run_diffsed_vs_shear
% Steps:
% 1. build both kernel fields
% 2. compare the weighted maps
% 3. save a simple summary

clear;
clc;

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(genpath(fullfile(repo_root, 'src')));

fig_dir = fullfile(repo_root, 'output', 'figures');
log_dir = fullfile(repo_root, 'output', 'logs');

if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end
if ~exist(log_dir, 'dir')
    mkdir(log_dir);
end

diam_cm = logspace(-4, 0, 220);
law_names = {'current', 'kriest_8', 'kriest_9', 'siegel_2025'};
eps_vals_mks = [1e-8, 1e-6, 1e-4];
conc_amp = 1e3;
conc_exp = -2.5;
mode_tag = 'direct_sinking_law';

out = compare_diffsed_vs_shear(diam_cm, law_names, eps_vals_mks, ...
    conc_amp, conc_exp, fig_dir, log_dir, mode_tag);

disp('Saved summary text:');
disp(out.summary_text);
disp('Summary rows:');
disp(out.summary);

function out = compare_diffsed_vs_shear(diam_cm, law_names, eps_vals_mks, conc_amp, conc_exp, fig_dir, log_dir, mode_tag)
diam_um = diam_cm .* 1e4;
[D1, D2] = ndgrid(diam_cm, diam_cm);

C = local_powerlaw_concentration(diam_cm, conc_amp, conc_exp);
[C1, C2] = ndgrid(C, C);

rows = struct('epsilon_mks', {}, 'law', {}, 'frac_ds_gt_shear', {}, ...
    'ds_max', {}, 'shear_max', {}, 'ratio_max', {}, ...
    'peak_ds_d1_um', {}, 'peak_ds_d2_um', {}, ...
    'peak_ratio_d1_um', {}, 'peak_ratio_d2_um', {});

for ie = 1:numel(eps_vals_mks)
    eps_now = eps_vals_mks(ie);
    tag = local_eps_tag(eps_now);

    ds_maps = cell(numel(law_names), 1);
    sh_maps = cell(numel(law_names), 1);
    ratio_maps = cell(numel(law_names), 1);
    all_ds = [];
    all_sh = [];
    all_ratio = [];

    for il = 1:numel(law_names)
        law = law_names{il};
        w = local_sinking_speed_named(diam_cm, law);
        [W1, W2] = ndgrid(w, w);

        beta_ds = local_beta_diff_sed(D1, D2, W1, W2);
        beta_sh = local_beta_turb_shear(D1, D2, eps_now);

        ds_rate = beta_ds .* C1 .* C2;
        sh_rate = beta_sh .* C1 .* C2;
        ratio = ds_rate ./ max(sh_rate, realmin);

        z_ds = log10(max(ds_rate, realmin));
        z_sh = log10(max(sh_rate, realmin));
        z_ratio = log10(max(ratio, realmin));

        ds_maps{il} = z_ds;
        sh_maps{il} = z_sh;
        ratio_maps{il} = z_ratio;

        all_ds = [all_ds; z_ds(isfinite(z_ds))]; %#ok<AGROW>
        all_sh = [all_sh; z_sh(isfinite(z_sh))]; %#ok<AGROW>
        all_ratio = [all_ratio; z_ratio(isfinite(z_ratio))]; %#ok<AGROW>

        [ds_max, ds_idx] = max(ds_rate(:));
        [ratio_max, ratio_idx] = max(ratio(:));
        shear_max = max(sh_rate(:));
        frac_ds = mean(ds_rate(:) > sh_rate(:));

        [id1_ds, id2_ds] = ind2sub(size(ds_rate), ds_idx);
        [id1_rt, id2_rt] = ind2sub(size(ratio), ratio_idx);

        row = struct();
        row.epsilon_mks = eps_now;
        row.law = string(law);
        row.frac_ds_gt_shear = frac_ds;
        row.ds_max = ds_max;
        row.shear_max = shear_max;
        row.ratio_max = ratio_max;
        row.peak_ds_d1_um = diam_um(id1_ds);
        row.peak_ds_d2_um = diam_um(id2_ds);
        row.peak_ratio_d1_um = diam_um(id1_rt);
        row.peak_ratio_d2_um = diam_um(id2_rt);
        rows(end + 1) = row; %#ok<AGROW>
    end

    ds_lim = [prctile(all_ds, 1), prctile(all_ds, 99)];
    sh_lim = [prctile(all_sh, 1), prctile(all_sh, 99)];
    rt_lim = [prctile(all_ratio, 1), prctile(all_ratio, 99)];

    local_plot_map_set(diam_um, law_names, ds_maps, ds_lim, ...
        sprintf('Differential settling (eps = %.1e)', eps_now), ...
        'log10(beta ds * C1 * C2)', fullfile(fig_dir, ['diffsed_c1c2_eps_' tag '.png']));

    local_plot_map_set(diam_um, law_names, sh_maps, sh_lim, ...
        sprintf('Turbulent shear (eps = %.1e)', eps_now), ...
        'log10(beta shear * C1 * C2)', fullfile(fig_dir, ['shear_c1c2_eps_' tag '.png']));

    local_plot_map_set(diam_um, law_names, ratio_maps, rt_lim, ...
        sprintf('DS / shear ratio (eps = %.1e, %s)', eps_now, mode_tag), ...
        'log10((beta ds * C1 * C2)/(beta shear * C1 * C2))', ...
        fullfile(fig_dir, ['diffsed_vs_shear_ratio_eps_' tag '.png']));
end

S = struct2table(rows);
summary_text = fullfile(log_dir, 'diffsed_vs_shear_summary.txt');
fid = fopen(summary_text, 'w');
fprintf(fid, 'Differential settling vs turbulent shear\n\n');
fprintf(fid, 'C(d) = %.1e * d^(%.1f)\n', conc_amp, conc_exp);
fprintf(fid, 'Size range: 1 micron to 1 cm\n\n');
fprintf(fid, 'Kernel mode used for these figures: %s\n', mode_tag);
fprintf(fid, 'Plotted ratio is DS/shear (not shear/DS).\n');
fprintf(fid, 'Color convention in figure: high value = DS dominant.\n\n');

for i = 1:height(S)
    fprintf(fid, 'epsilon = %.1e | %s\n', S.epsilon_mks(i), char(S.law(i)));
    fprintf(fid, '  frac(ds > shear) = %.3f\n', S.frac_ds_gt_shear(i));
    fprintf(fid, '  ds max = %.3e at (%.0f um, %.0f um)\n', ...
        S.ds_max(i), S.peak_ds_d1_um(i), S.peak_ds_d2_um(i));
    fprintf(fid, '  shear max = %.3e\n', S.shear_max(i));
    fprintf(fid, '  ratio max = %.3e at (%.0f um, %.0f um)\n\n', ...
        S.ratio_max(i), S.peak_ratio_d1_um(i), S.peak_ratio_d2_um(i));
end
fclose(fid);

out = struct();
out.summary = S;
out.summary_text = summary_text;
end

function local_plot_map_set(diam_um, law_names, maps, clim, fig_title, cbar_label, out_path)
labels = {'current', 'kriest 8', 'kriest 9', 'siegel 2025'};
fig = figure('Color', 'w', 'Position', [110 110 900 730]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:numel(law_names)
    ax = nexttile(tl, i);
    z = maps{i};
    z = min(max(z, clim(1)), clim(2));
    imagesc(ax, log10(diam_um), log10(diam_um), z);
    axis(ax, 'xy');
    caxis(ax, clim);
    xlabel(ax, 'log10 d1 (um)');
    ylabel(ax, 'log10 d2 (um)');
    title(ax, labels{i});
    ax.LineWidth = 1.0;
    ax.FontSize = 11;
end
colormap(fig, parula);
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = cbar_label;
title(tl, fig_title);
save_figure(fig, out_path);
end

function tag = local_eps_tag(eps_now)
tag = strrep(sprintf('%.0e', eps_now), 'e-0', 'em0');
tag = strrep(tag, 'e-', 'em');
end

function c = local_powerlaw_concentration(d_cm, amp, expo)
ref = min(d_cm(d_cm > 0));
if isempty(ref)
    ref = 1.0;
end
c = amp .* (d_cm ./ ref) .^ expo;
c(~isfinite(c)) = 0;
c(c < 0) = 0;
end

function beta = local_beta_diff_sed(d1_cm, d2_cm, w1, w2)
beta = (pi / 4.0) .* (d1_cm + d2_cm) .* (d1_cm + d2_cm) .* abs(w1 - w2);
end

function beta = local_beta_turb_shear(d1_cm, d2_cm, eps_mks)
cfg = SimulationConfig();
r1 = d1_cm ./ (2.0 * cfg.r_to_rg);
r2 = d2_cm ./ (2.0 * cfg.r_to_rg);
p = min(r1, r2) ./ max(r1, r2);
p1 = 1.0 + p;
eff = 1.0 - (1.0 + 5.0 .* p + 2.5 .* p .* p) ./ (p1 .^ 5);
rg = (r1 + r2) .* cfg.r_to_rg;
shape = sqrt(8.0 * pi / 15.0) .* eff .* rg .^ 3;
eps_cgs = eps_mks * 1e4;
gamma = sqrt(eps_cgs / cfg.kvisc);
beta = shape .* gamma;
beta(~isfinite(beta)) = 0;
beta(beta < 0) = 0;
end

function w_cm_s = local_sinking_speed_named(d_cm, law_name)
cfg = SimulationConfig();
law = lower(string(law_name));
switch law
    case "current"
        w_cm_s = local_current_law(d_cm, cfg);
    case "kriest_8"
        w_cm_s = (66 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "kriest_9"
        w_cm_s = (132 .* (d_cm .^ 0.62) .* 100) ./ cfg.day_to_sec;
    case "siegel_2025"
        d_mm = d_cm .* 10.0;
        w_cm_s = (20.2 .* (d_mm .^ 0.67) .* 100) ./ cfg.day_to_sec;
    otherwise
        error('Unknown law: %s', char(law_name));
end
w_cm_s(~isfinite(w_cm_s)) = 0;
w_cm_s(w_cm_s < 0) = 0;
end

function w_cm_s = local_current_law(d_cm, cfg)
r_v = 0.5 .* d_cm;
setcon = KernelLibrary.currentSetcon(cfg);
r_i = KernelLibrary.conservativeToFractalRadius(r_v, cfg);
w_cm_s = setcon .* (r_v .^ 3) ./ max(r_i, realmin);
end

function save_figure(fig_handle, fig_path)
[fig_parent, ~, ~] = fileparts(fig_path);
if ~exist(fig_parent, 'dir')
    mkdir(fig_parent);
end
set(fig_handle, 'PaperPositionMode', 'auto');
try
    exportgraphics(fig_handle, fig_path, 'Resolution', 220);
catch
    saveas(fig_handle, fig_path);
end
end

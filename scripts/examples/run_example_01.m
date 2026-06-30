% run_example_01.m
% Example 1: pulse at top, constant eps, sinking + coagulation only.
%
% Put a power-law size distribution into the top layer at t=0.
% Run 30 days with constant turbulence. No zoo, no disagg.
% Check: particles sink, large bins gain mass, no negatives.
%
% Saves: docs/figures/example_01_profile.png

script_dir = fileparts(mfilename('fullpath'));
addpath(script_dir);
addpath(fullfile(script_dir, '..', '..', 'src'));

set(0,'DefaultAxesFontName','Arial');
set(0,'DefaultTextFontName','Arial');

% --- grid ---
col_grid  = ColumnGrid(1000, 20);
z_centers = col_grid.z_centers;

% --- config ---
cfg = SimulationConfig();
cfg.n_sections     = 30;
cfg.sinking_law    = 'kriest_8';
cfg.ds_kernel_mode = 'sinking_law';
cfg.r_to_rg        = 1.6;
cfg.alpha          = 1.0;
cfg.enable_coag    = true;
cfg.enable_disagg  = false;
cfg.enable_zoo     = false;
cfg.enable_microbe = false;

% constant turbulence profile: use typical() then override eps
prof = DepthProfile.typical(col_grid.z_centers);
prof.eps(:) = 1e-5;   % cm^2/s^3, uniform mid-water value

sim  = ColumnSimulation(cfg, col_grid, prof);

% power-law IC in top layer: phi ~ d^(-4), normalized to BV = 1e-6
d_cm   = sim.size_grid.dcomb(:)';
phi_pl = d_cm .^ (-4);
phi_pl = phi_pl / sum(phi_pl) * 1e-6;

Y   = zeros(col_grid.n_z, cfg.n_sections);
Yfp = zeros(col_grid.n_z, cfg.n_sections);
Y(1,:) = phi_pl;

% save profiles at day 0, 10, 30
save_days = [0, 10, 30];
profiles  = zeros(col_grid.n_z, numel(save_days));
profiles(:,1) = sum(Y, 2);

dt    = 0.25;
n_day = 30;
steps = round(n_day / dt);
day_now = 0;

for k = 1:steps
    [Y, Yfp] = sim.rhs.stepY(Y, dt, Yfp);
    day_now = day_now + dt;
    for s = 2:numel(save_days)
        if abs(day_now - save_days(s)) < dt/2
            profiles(:,s) = sum(Y + Yfp, 2);
        end
    end
end

% check no negatives
if any(Y(:) < 0) || any(Yfp(:) < 0)
    warning('Negative concentrations found!');
else
    fprintf('No negatives. Good.\n');
end

% --- figure ---
fs = 8;
ls = {'k:', 'k--', 'k-'};
lw = [1.0, 1.2, 1.4];

figure('Units','centimeters','Position',[2 2 14 8],'Color','white');

% left: total BV at 3 times (log x)
subplot(1,2,1);
hold on;
for s = 1:numel(save_days)
    pv = profiles(:,s);
    pv(pv <= 0) = NaN;
    plot(pv, z_centers, ls{s}, 'LineWidth', lw(s), ...
        'DisplayName', sprintf('t=%dd', save_days(s)));
end
set(gca, 'YDir','reverse', 'XScale','log', 'FontSize',fs, 'Box','on', ...
    'YLim',[0 1000]);
xlabel('BV (m^3 m^{-3})', 'FontSize',fs);
ylabel('depth (m)', 'FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('(a) total', 'FontWeight','normal','FontSize',fs);

% right: by size group at t=30 (log x)
subplot(1,2,2);
hold on;
bin_groups = {1:5, 10:15, 20:25};
grp_labels = {'20-50 \mum', '0.3-1 mm', '4-16 mm'};
for g = 1:3
    pv = sum(Y(:, bin_groups{g}), 2);
    pv(pv <= 0) = NaN;
    plot(pv, z_centers, ls{g}, 'LineWidth', lw(g), ...
        'DisplayName', grp_labels{g});
end
set(gca, 'YDir','reverse', 'XScale','log', 'FontSize',fs, 'Box','on', ...
    'YLim',[0 1000], 'YTickLabel',{});
xlabel('BV (m^3 m^{-3})', 'FontSize',fs);
legend('Location','southeast','FontSize',fs,'Box','off');
title('(b) by size, t=30d', 'FontWeight','normal','FontSize',fs);

fig_dir = fullfile(script_dir, '..', '..', 'docs', 'figures');
if ~exist(fig_dir,'dir'), mkdir(fig_dir); end
print(gcf, fullfile(fig_dir,'example_01_profile.png'), '-dpng', '-r150');
fprintf('Saved example_01_profile.png\n');

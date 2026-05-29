function out = validate_coagulation_conservation(sim_base, sim_coag, fig_dir, log_dir, tab_dir, tag)
% validate_coagulation_conservation
% Make simple trust plots and summary for step-4.

t_day = sim_base.t_s(:) ./ 86400.0;
size_um = sim_base.size_um(:);

v0 = max(sim_base.tracked_volume_total(1), realmin);
err_base = 100.0 .* (sim_base.tracked_volume_total - sim_base.tracked_volume_total(1)) ./ ...
    max(abs(sim_base.tracked_volume_total(1)), realmin);
err_coag = 100.0 .* (sim_coag.tracked_volume_total - sim_coag.tracked_volume_total(1)) ./ ...
    max(abs(sim_coag.tracked_volume_total(1)), realmin);

small_mask = size_um < 500;
large_mask = size_um >= 500;
base_large = 100.0 .* sum(sim_base.column_volume_by_size(:, large_mask), 2) ./ v0;
coag_large = 100.0 .* sum(sim_coag.column_volume_by_size(:, large_mask), 2) ./ v0;

% conservation
fig1 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax1 = axes(fig1); hold(ax1, 'on');
plot(ax1, t_day, err_base, 'k', 'LineWidth', 1.4, 'DisplayName', 'no coag');
plot(ax1, t_day, err_coag, 'r', 'LineWidth', 1.4, 'DisplayName', 'with coag');
xlabel(ax1, 'Time (day)');
ylabel(ax1, 'Tracked volume error (%)');
title(ax1, 'Coagulation conservation');
legend(ax1, 'Location', 'best', 'Box', 'off');
ax1.LineWidth = 1.0; ax1.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, [tag '_conservation.png']));
close(fig1);

% large size volume
fig2 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax2 = axes(fig2); hold(ax2, 'on');
plot(ax2, t_day, base_large, 'k', 'LineWidth', 1.4, 'DisplayName', 'no coag');
plot(ax2, t_day, coag_large, 'r', 'LineWidth', 1.4, 'DisplayName', 'with coag');
xlabel(ax2, 'Time (day)');
ylabel(ax2, 'Volume in sizes >= 500 um (%)');
title(ax2, 'Large-size volume');
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0; ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, [tag '_large_size_volume.png']));
close(fig2);

% column final PSD
fig3 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax3 = axes(fig3); hold(ax3, 'on');
plot(ax3, size_um, max(sim_base.column_number(end, :)', realmin), 'k', 'LineWidth', 1.4, 'DisplayName', 'no coag');
plot(ax3, size_um, max(sim_coag.column_number(end, :)', realmin), 'r', 'LineWidth', 1.4, 'DisplayName', 'with coag');
set(ax3, 'XScale', 'log', 'YScale', 'log');
xlabel(ax3, 'Particle size (um)');
ylabel(ax3, 'Particles left in column');
title(ax3, 'Final column PSD');
legend(ax3, 'Location', 'best', 'Box', 'off');
ax3.LineWidth = 1.0; ax3.FontSize = 11;
save_figure(fig3, fullfile(fig_dir, [tag '_column_psd.png']));
close(fig3);

% final depth PSD (all sizes)
fig4 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax4 = axes(fig4);
imagesc(ax4, size_um, sim_base.z_m, sim_coag.conc(:, :, end));
axis(ax4, 'xy');
set(ax4, 'YDir', 'reverse', 'XScale', 'log');
xlabel(ax4, 'Particle size (um)');
ylabel(ax4, 'Depth (m)');
title(ax4, 'Final depth-size field (with coag)');
colormap(fig4, parula);
cb = colorbar(ax4);
cb.Label.String = 'Concentration';
ax4.LineWidth = 1.0; ax4.FontSize = 11;
save_figure(fig4, fullfile(fig_dir, [tag '_final_depth_psd.png']));
close(fig4);

% bottom signal through time
fig5 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax5 = axes(fig5);
imagesc(ax5, t_day, size_um, sim_coag.bottom_signal');
axis(ax5, 'xy');
set(ax5, 'YScale', 'log');
xlabel(ax5, 'Time (day)');
ylabel(ax5, 'Particle size (um)');
title(ax5, 'Bottom signal by size');
colormap(fig5, parula);
cb = colorbar(ax5);
cb.Label.String = 'Bottom concentration';
ax5.LineWidth = 1.0; ax5.FontSize = 11;
save_figure(fig5, fullfile(fig_dir, [tag '_bottom_psd_time.png']));
close(fig5);

% table
summary = table(size_um, ...
    sim_base.column_number(end, :)', sim_coag.column_number(end, :)', ...
    sim_coag.column_number(end, :)' ./ max(sim_base.column_number(end, :)', realmin), ...
    'VariableNames', {'size_um', 'base_final_column_number', 'coag_final_column_number', 'coag_to_base_ratio'});
csv_path = fullfile(tab_dir, [tag '_summary.csv']);
writetable(summary, csv_path);

% log
log_path = fullfile(log_dir, [tag '.txt']);
fid = fopen(log_path, 'w');
fprintf(fid, 'Coagulation validation\n\n');
fprintf(fid, 'Main checks:\n');
fprintf(fid, '- neg_count (no coag) = %d\n', sum(sim_base.conc(:) < -1e-12));
fprintf(fid, '- neg_count (with coag) = %d\n', sum(sim_coag.conc(:) < -1e-12));
fprintf(fid, '- max tracked-volume error (no coag) = %.6e %%\n', max(abs(err_base)));
fprintf(fid, '- max tracked-volume error (with coag) = %.6e %%\n', max(abs(err_coag)));
fprintf(fid, '- final total-number change (no coag) = %.6f %%\n', ...
    100.0 .* (sim_base.total_number(end) - sim_base.total_number(1)) ./ max(abs(sim_base.total_number(1)), realmin));
fprintf(fid, '- final total-number change (with coag) = %.6f %%\n', ...
    100.0 .* (sim_coag.total_number(end) - sim_coag.total_number(1)) ./ max(abs(sim_coag.total_number(1)), realmin));
fclose(fid);

out = struct('csv_path', csv_path, 'log_path', log_path);
end


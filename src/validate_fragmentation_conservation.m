function out = validate_fragmentation_conservation(sim_coag, sim_frag, fig_dir, log_dir, tab_dir, tag)
% validate_fragmentation_conservation
% Make simple trust plots and summary for step-5.

t_day = sim_coag.t_s(:) ./ 86400.0;
size_um = sim_coag.size_um(:);

err_coag = 100.0 .* (sim_coag.tracked_volume_total - sim_coag.tracked_volume_total(1)) ./ ...
    max(abs(sim_coag.tracked_volume_total(1)), realmin);
err_frag = 100.0 .* (sim_frag.tracked_volume_total - sim_frag.tracked_volume_total(1)) ./ ...
    max(abs(sim_frag.tracked_volume_total(1)), realmin);

small_mask = size_um <= 500;
v0 = max(sim_coag.tracked_volume_total(1), realmin);
coag_small = 100.0 .* sum(sim_coag.column_volume_by_size(:, small_mask), 2) ./ v0;
frag_small = 100.0 .* sum(sim_frag.column_volume_by_size(:, small_mask), 2) ./ v0;

% conservation
fig1 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax1 = axes(fig1); hold(ax1, 'on');
plot(ax1, t_day, err_coag, 'k', 'LineWidth', 1.4, 'DisplayName', 'coag only');
plot(ax1, t_day, err_frag, 'r', 'LineWidth', 1.4, 'DisplayName', 'coag + frag');
xlabel(ax1, 'Time (day)');
ylabel(ax1, 'Tracked volume error (%)');
title(ax1, 'Fragmentation conservation');
legend(ax1, 'Location', 'best', 'Box', 'off');
ax1.LineWidth = 1.0; ax1.FontSize = 11;
save_figure(fig1, fullfile(fig_dir, [tag '_conservation.png']));
close(fig1);

% small size volume
fig2 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax2 = axes(fig2); hold(ax2, 'on');
plot(ax2, t_day, coag_small, 'k', 'LineWidth', 1.4, 'DisplayName', 'coag only');
plot(ax2, t_day, frag_small, 'r', 'LineWidth', 1.4, 'DisplayName', 'coag + frag');
xlabel(ax2, 'Time (day)');
ylabel(ax2, 'Volume in sizes <= 500 um (%)');
title(ax2, 'Small-size volume');
legend(ax2, 'Location', 'best', 'Box', 'off');
ax2.LineWidth = 1.0; ax2.FontSize = 11;
save_figure(fig2, fullfile(fig_dir, [tag '_small_size_volume.png']));
close(fig2);

% column final PSD
fig3 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax3 = axes(fig3); hold(ax3, 'on');
plot(ax3, size_um, max(sim_coag.column_number(end, :)', realmin), 'k', 'LineWidth', 1.4, 'DisplayName', 'coag only');
plot(ax3, size_um, max(sim_frag.column_number(end, :)', realmin), 'r', 'LineWidth', 1.4, 'DisplayName', 'coag + frag');
set(ax3, 'XScale', 'log', 'YScale', 'log');
xlabel(ax3, 'Particle size (um)');
ylabel(ax3, 'Particles left in column');
title(ax3, 'Final column PSD');
legend(ax3, 'Location', 'best', 'Box', 'off');
ax3.LineWidth = 1.0; ax3.FontSize = 11;
save_figure(fig3, fullfile(fig_dir, [tag '_column_psd.png']));
close(fig3);

% final depth PSD
fig4 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax4 = axes(fig4);
imagesc(ax4, size_um, sim_coag.z_m, sim_frag.conc(:, :, end));
axis(ax4, 'xy');
set(ax4, 'YDir', 'reverse', 'XScale', 'log');
xlabel(ax4, 'Particle size (um)');
ylabel(ax4, 'Depth (m)');
title(ax4, 'Final depth-size field (coag + frag)');
colormap(fig4, parula);
cb = colorbar(ax4);
cb.Label.String = 'Concentration';
ax4.LineWidth = 1.0; ax4.FontSize = 11;
save_figure(fig4, fullfile(fig_dir, [tag '_final_depth_psd.png']));
close(fig4);

% bottom signal through time
fig5 = figure('Color', 'w', 'Position', [120 120 700 520]);
ax5 = axes(fig5);
imagesc(ax5, t_day, size_um, sim_frag.bottom_signal');
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
    sim_coag.column_number(end, :)', sim_frag.column_number(end, :)', ...
    sim_frag.column_number(end, :)' ./ max(sim_coag.column_number(end, :)', realmin), ...
    'VariableNames', {'size_um', 'coag_final_column_number', 'frag_final_column_number', 'frag_to_coag_ratio'});
csv_path = fullfile(tab_dir, [tag '_summary.csv']);
writetable(summary, csv_path);

% log
log_path = fullfile(log_dir, [tag '.txt']);
fid = fopen(log_path, 'w');
fprintf(fid, 'Fragmentation validation\n\n');
fprintf(fid, 'Main checks:\n');
fprintf(fid, '- neg_count (coag only) = %d\n', sum(sim_coag.conc(:) < -1e-12));
fprintf(fid, '- neg_count (coag + frag) = %d\n', sum(sim_frag.conc(:) < -1e-12));
fprintf(fid, '- max tracked-volume error (coag only) = %.6e %%\n', max(abs(err_coag)));
fprintf(fid, '- max tracked-volume error (coag + frag) = %.6e %%\n', max(abs(err_frag)));
fprintf(fid, '- final total-number change (coag only) = %.6f %%\n', ...
    100.0 .* (sim_coag.total_number(end) - sim_coag.total_number(1)) ./ max(abs(sim_coag.total_number(1)), realmin));
fprintf(fid, '- final total-number change (coag + frag) = %.6f %%\n', ...
    100.0 .* (sim_frag.total_number(end) - sim_frag.total_number(1)) ./ max(abs(sim_frag.total_number(1)), realmin));
fclose(fid);

out = struct('csv_path', csv_path, 'log_path', log_path);
end


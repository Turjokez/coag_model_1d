function save_figure(fig_handle, fig_path)
% save_figure
% Short note:
% 1. make sure folder exists
% 2. save figure with stable settings

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


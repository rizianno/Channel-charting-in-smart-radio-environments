function plot_aoa_error_map(geo, err_deg, outFile, cmax, keepOpen)
% PLOT_AOA_ERROR_MAP  Map of the per-user joint AoA error, with EMS and BS.
%   geo     : struct with R_xy [U x 2] (evaluated users), EMS_position [2 x 3],
%             BS_position [1 x 3]
%   err_deg : [U x 1] joint AoA error of the same users [deg]
%   outFile : PNG path
%   cmax    : top of the color scale [deg]; larger errors take the top color
%   keepOpen: true = leave the figure on screen (interactive use)
    if nargin < 4, cmax = 30; end
    if nargin < 5, keepOpen = false; end
    f = figure('Position',[100 100 700 600]); hold on; grid on; box on; axis equal;
    if keepOpen, f.Visible = 'on'; end
    scatter(geo.R_xy(:,1), geo.R_xy(:,2), 12, min(err_deg, cmax), 'filled');
    clim([0 cmax]);
    cb_ = colorbar; ylabel(cb_, 'Joint AoA error [deg]');
    steps = [1 2 5 10 20 50 100];
    step  = steps(find(cmax./steps <= 6, 1));
    cb_.Ticks = unique([0:step:cmax, cmax]);
    cb_.TickLabels = compose('%g', cb_.Ticks);
    cb_.TickLabels{end} = sprintf('\\geq%g', cmax);
    plot(geo.EMS_position(:,1), geo.EMS_position(:,2), 'ks', 'MarkerSize',10, 'MarkerFaceColor','y');
    plot(geo.BS_position(1), geo.BS_position(2), 'k^', 'MarkerSize',10, 'MarkerFaceColor','r');
    text(geo.EMS_position(:,1), geo.EMS_position(:,2)-3, {'EMS1','EMS2'}, 'HorizontalAlignment','center');
    text(geo.BS_position(1), geo.BS_position(2)+3, 'BS', 'HorizontalAlignment','center');
    xlabel('x [m]'); ylabel('y [m]');
    title(sprintf('Spatial distribution of the AoA error (local fusion), %d users', numel(err_deg)));
    exportgraphics(f, outFile, 'Resolution',150);
    if ~keepOpen, close(f); end
end

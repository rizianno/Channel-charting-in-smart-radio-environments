function plot_aoa_error_map(geo, err_deg, outFile, cmax)
% PLOT_AOA_ERROR_MAP  Map of the per-user joint AoA error, with EMS and BS.
%   geo     : struct with R_xy [U x 2] (evaluated users), EMS_position [2 x 3],
%             BS_position [1 x 3]
%   err_deg : [U x 1] joint AoA error of the same users [deg]
%   outFile : PNG path
%   cmax    : top of the color scale [deg]; larger errors take the top color
    if nargin < 4, cmax = 40; end
    f = figure('Position',[100 100 700 600]); hold on; grid on; box on; axis equal;
    scatter(geo.R_xy(:,1), geo.R_xy(:,2), 12, min(err_deg, cmax), 'filled');
    clim([0 cmax]);
    cb_ = colorbar; ylabel(cb_, 'Joint AoA error [deg]');
    cb_.Ticks = 0:10:cmax;
    cb_.TickLabels{end} = sprintf('\\geq%g', cmax);
    plot(geo.EMS_position(:,1), geo.EMS_position(:,2), 'ks', 'MarkerSize',10, 'MarkerFaceColor','y');
    plot(geo.BS_position(1), geo.BS_position(2), 'k^', 'MarkerSize',10, 'MarkerFaceColor','r');
    text(geo.EMS_position(:,1), geo.EMS_position(:,2)-3, {'EMS1','EMS2'}, 'HorizontalAlignment','center');
    text(geo.BS_position(1), geo.BS_position(2)+3, 'BS', 'HorizontalAlignment','center');
    xlabel('x [m]'); ylabel('y [m]');
    title(sprintf('Spatial distribution of the AoA error (local fusion), %d users', numel(err_deg)));
    exportgraphics(f, outFile, 'Resolution',150); close(f);
end

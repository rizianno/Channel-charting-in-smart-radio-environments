function plot_aoa_error_map_per_ems(geo, Y_true, Y_est, outFile, cmax, mode, keepOpen)
% PLOT_AOA_ERROR_MAP_PER_EMS  Two maps side by side: AoA error seen from EMS1
% and from EMS2, each with its own colorbar and range.
%   geo      : struct with R_xy [U x 2], EMS_position [2 x >=2], BS_position [1 x >=2]
%   Y_true   : [U x 2] true azimuth from EMS1, EMS2 [deg]
%   Y_est    : [U x 2] estimated azimuth [deg]
%   outFile  : PNG path
%   cmax     : [1 x 2] (or scalar) top of the color scale per EMS [deg]
%   mode     : 'abs'    |est - true|, 0..cmax (default)
%              'signed' est - true, -cmax..cmax (blue = under-, red = over-estimate)
%   keepOpen : true = leave the figure on screen (interactive use)
    if nargin < 5 || isempty(cmax), cmax = 30; end
    if nargin < 6 || isempty(mode), mode = 'abs'; end
    if nargin < 7, keepOpen = false; end
    if isscalar(cmax), cmax = [cmax cmax]; end
    signed = strcmpi(mode, 'signed');

    f = figure('Position',[100 100 1300 600]);
    if keepOpen, f.Visible = 'on'; end
    tl = tiledlayout(f, 1, 2, 'TileSpacing','compact', 'Padding','compact');

    for e = 1:2
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
        d = Y_est(:,e) - Y_true(:,e);
        c = cmax(e);
        if signed
            v = max(min(d, c), -c);
            clim(ax, [-c c]); colormap(ax, diverging_map(256));
            lbl = sprintf('EMS%d AoA error, estimate - true [deg]', e);
        else
            v = min(abs(d), c);
            clim(ax, [0 c]); colormap(ax, parula(256));
            lbl = sprintf('EMS%d AoA error |estimate - true| [deg]', e);
        end
        scatter(ax, geo.R_xy(:,1), geo.R_xy(:,2), 12, v, 'filled');

        cb_ = colorbar(ax); ylabel(cb_, lbl);
        [cb_.Ticks, cb_.TickLabels] = ticks_for(c, signed);

        % both surfaces and the BS; the one this map refers to is highlighted
        other = 3 - e;
        plot(ax, geo.EMS_position(other,1), geo.EMS_position(other,2), 'ks', ...
             'MarkerSize',9, 'MarkerFaceColor',[0.8 0.8 0.8]);
        plot(ax, geo.EMS_position(e,1), geo.EMS_position(e,2), 'ks', ...
             'MarkerSize',12, 'LineWidth',1.5, 'MarkerFaceColor','y');
        plot(ax, geo.BS_position(1), geo.BS_position(2), 'k^', 'MarkerSize',10, 'MarkerFaceColor','r');
        text(ax, geo.EMS_position(:,1), geo.EMS_position(:,2)-3, {'EMS1','EMS2'}, ...
             'HorizontalAlignment','center');
        text(ax, geo.BS_position(1), geo.BS_position(2)+3, 'BS', 'HorizontalAlignment','center');
        xlabel(ax,'x [m]'); ylabel(ax,'y [m]');
        title(ax, sprintf('Seen from EMS%d: median %.1f%s, %.0f%% of users within %g%s', ...
              e, median(abs(d)), char(176), 100*mean(abs(d) < c), c, char(176)));
    end
    if signed, kind = 'signed'; else, kind = 'absolute'; end
    title(tl, sprintf('AoA error per EMS (%s, local fusion), %d users', kind, size(Y_true,1)));

    exportgraphics(f, outFile, 'Resolution',150);
    if ~keepOpen, close(f); end
end

function [t, lbl] = ticks_for(c, signed)
% readable ticks up to c; the end label(s) mark saturation (>= c)
    steps = [1 2 5 10 20 50 100];
    step  = steps(find(c./steps <= 6, 1));
    t = unique([0:step:c, c]);
    if signed
        t = unique([-fliplr(t), t]) + 0;   % + 0 turns -0 into 0
    end
    lbl = compose('%g', t);
    lbl{end} = sprintf('\\geq%g', c);
    if signed, lbl{1} = sprintf('\\leq-%g', c); end
end

function m = diverging_map(n)
% blue - white - red
    a = [0.23 0.30 0.75]; w = [0.97 0.97 0.97]; b = [0.71 0.02 0.15];
    h = floor(n/2); s = linspace(0,1,h)'; s2 = linspace(0,1,n-h)';
    m = [a + s.*(w-a); w + s2.*(b-w)];
end

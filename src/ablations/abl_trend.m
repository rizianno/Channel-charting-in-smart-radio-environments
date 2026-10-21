function abl_trend(T, scn, outDir)
% ABL_TREND  Error vs codebook size of one scenario (all K): trend.csv and trend.png.
% Plot: median |error| on the natural subset vs number of acquisitions, per method;
% for F0 the gap between typical_single and the fused lines is the fusion gain.
    T = T(~startsWith(T.method,'single_'), :);
    T.method = regexprep(T.method, '^best_single\(cfg\d+\)$', 'best_single');
    writetable(T, fullfile(outDir,'trend.csv'));

    panels = {'EMS1','EMS2','pooled'};
    f = figure('Visible','off', 'Position',[100 100 1300 400]);
    tl = tiledlayout(f, 1, 3, 'TileSpacing','compact');
    meths = unique(T.method, 'stable');
    for p = 1:3
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for m = meths.'
            r = T(strcmp(T.method,m{1}) & strcmp(T.subset,'natural') & strcmp(T.ems,panels{p}), :);
            r = sortrows(r, 'acquisitions');
            ls = '-o'; if contains(m{1},'single'), ls = ':s'; end
            plot(ax, r.acquisitions, r.median_deg, ls, 'LineWidth',1.4, 'DisplayName', strrep(m{1},'_',' '));
        end
        set(ax, 'XScale','log'); xlabel(ax, 'acquisitions (EMS configurations measured)');
        ylabel(ax, 'median |AoA error| [deg]'); title(ax, panels{p});
        if p == 1, legend(ax, 'Location','northeast'); end
    end
    title(tl, strrep(sprintf('%s: error vs codebook size (natural subset)', scn.name), '_', '\_'));
    exportgraphics(f, fullfile(outDir,'trend.png'), 'Resolution', 110); close(f);
end

function abl_compare(T, outDir, groupBy)
% ABL_COMPARE  Cross-scenario comparison of the main method of each scenario
% (F0: local_orig and local_free; F1-F3: chart): comparison.csv and one figure per subset.
%   groupBy : 'structure' (one panel per structure) or 'noise' (calibration)
    keep = ismember(T.method, {'local_orig','local_free','chart'});
    T = T(keep,:);
    writetable(T, fullfile(outDir,'comparison.csv'));
    for s = {'common','natural'}
        f = figure('Visible','off', 'Position',[100 100 1500 460]);
        if strcmp(groupBy,'noise'), G = T.noise; else, G = extractBefore(T.scenario,'_'); end
        groups = unique(G, 'stable');
        tl = tiledlayout(f, 1, numel(groups), 'TileSpacing','compact');
        for g = groups.'
            ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
            r = T(strcmp(G,g{1}) & strcmp(T.subset,s{1}) & strcmp(T.ems,'pooled'), :);
            lines = unique(strcat(r.scenario, '|', r.method), 'stable');
            for L = lines.'
                q = sortrows(r(strcmp(strcat(r.scenario,'|',r.method), L{1}), :), 'acquisitions');
                lbl = strrep(regexprep(L{1}, '^[^_]+_', ''), '_', ' ');
                plot(ax, q.acquisitions, q.median_deg, '-o', 'LineWidth',1.3, 'DisplayName', lbl);
            end
            set(ax, 'XScale','log'); xlabel(ax, 'acquisitions'); ylabel(ax, 'pooled median |error| [deg]');
            title(ax, strrep(g{1},'_','\_')); legend(ax, 'Location','northeast', 'FontSize',7);
        end
        title(tl, sprintf('Pooled EMS1+EMS2 median AoA error, subset "%s"', s{1}));
        exportgraphics(f, fullfile(outDir, sprintf('comparison_%s.png', s{1})), 'Resolution', 110);
        close(f);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AoA ERROR vs CODEBOOK SIZE: median, mean, 90th percentile and RMSE for every ablation
% scenario of one noise model, from the per-user CSVs (no re-simulation).
%
% Users: the "common" set (in front of both surfaces, anchor in no structure), so every
% curve uses the same users. Errors: joint sqrt(e1^2+e2^2) per user (also for the
% per-surface charts, from their two 1-D estimates) and pooled EMS1+EMS2 per-angle errors.
%   ABL_NOISE  noise folder under results_aoa/ablations (default est_G60dB)
% Outputs in results_aoa/ablations/<noise>/: trend_all.csv, trend_<metric>.png
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root  = fileparts(fileparts(mfilename('fullpath')));
noise = getenv('ABL_NOISE');  if isempty(noise), noise = 'est_G60dB'; end
nDir  = fullfile(root,'results_aoa','ablations',noise);

d = dir(nDir);  scns = {d([d.isdir] & ~startsWith({d.name},'.')).name};
common = [];
rows = {};
for s = scns
    for K = 1:5
        f = dir(fullfile(nDir, s{1}, sprintf('K%d',K), 'per_user_*.csv'));
        if isempty(f), continue; end
        P = readtable(fullfile(f(1).folder, f(1).name));
        if isempty(common), common = common_users(P); end
        e1 = P.err1_deg(common);  e2 = P.err2_deg(common);
        sets = {'joint', sqrt(e1.^2 + e2.^2); 'pooled', abs([e1; e2])};
        for i = 1:2
            x = sets{i,2};
            rows(end+1,:) = {s{1}, K, K^2, sets{i,1}, median(x), mean(x), prctile(x,90), sqrt(mean(x.^2))}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'scenario','K','acquisitions','error','median_deg','mean_deg','p90_deg','rmse_deg'});
writetable(T, fullfile(nDir,'trend_all.csv'));
fprintf('%d scenarios, %d common users\n', numel(unique(T.scenario)), nnz(common));

% ---- figures: one per metric; rows = joint / pooled error; columns = user/chart structure
groups = {'S1-joint', 'joint chart, users in front of both surfaces'; ...
          'S2-perEMS', 'one chart per surface'; ...
          'S3-joint-allUE', 'joint chart, all users in t-SNE'};
feats  = {'F0-cov-per-sweep_dir-on',  'per-sweep charts + fusion, direct kept'; ...
          'F0-cov-per-sweep_dir-off', 'per-sweep charts + fusion, direct removed'; ...
          'F1-stacked-cov_dir-on',    'all sweeps stacked, direct kept'; ...
          'F1-stacked-cov_dir-off',   'all sweeps stacked, direct removed'; ...
          'F2-sweep-change_dir-on',   'sweep changes, direct kept'; ...
          'F2-sweep-change_dir-off',  'sweep changes, direct removed'; ...
          'F3-off-subtract_dir-on',   'surfaces-off subtracted, direct kept'};
cols = lines(size(feats,1));
metrics = {'median_deg','median'; 'mean_deg','mean'; 'p90_deg','90th percentile'; 'rmse_deg','RMSE'};
errs = {'joint','joint error \surd(e_1^2+e_2^2)'; 'pooled','per-surface error (EMS1 and EMS2 pooled)'};
for m = 1:size(metrics,1)
    f = figure('Position',[50 50 1700 900]);
    tl = tiledlayout(f, 2, size(groups,1), 'TileSpacing','compact');
    for r = 1:2
        for g = 1:size(groups,1)
            ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
            for k = 1:size(feats,1)
                for cb = {'', '_cb-wide'}
                    name = sprintf('%s_%s%s', groups{g,1}, feats{k,1}, cb{1});
                    q = sortrows(T(strcmp(T.scenario,name) & strcmp(T.error,errs{r,1}),:), 'acquisitions');
                    if isempty(q), continue; end
                    ls = '-o'; lbl = [feats{k,2} ', narrow'];
                    if ~isempty(cb{1}), ls = '--s'; lbl = [feats{k,2} ', wide']; end
                    plot(ax, q.acquisitions, q.(metrics{m,1}), ls, 'Color',cols(k,:), 'LineWidth',1.3, ...
                         'MarkerSize',4, 'DisplayName', lbl);
                end
            end
            set(ax,'XScale','log'); xticks(ax,[1 4 9 16 25]); xlim(ax,[0.9 28]);
            xlabel(ax,'EMS configurations (K^2)'); ylabel(ax, sprintf('%s [deg]', metrics{m,2}));
            title(ax, sprintf('%s: %s', groups{g,2}, errs{r,1}));
            if r == 1 && g == size(groups,1), legend(ax,'Location','eastoutside','FontSize',7); end
        end
    end
    title(tl, sprintf('Direct AoA estimation (%s): %s of the AoA error vs codebook size, %d common users', ...
          strrep(noise,'_','\_'), metrics{m,2}, nnz(common)));
    exportgraphics(f, fullfile(nDir, sprintf('trend_%s.png', strrep(metrics{m,2},' ','_'))), 'Resolution',110);
    close(f);
end
fprintf('=== DONE -> %s ===\n', nDir);

function c = common_users(P)
% in front of both surfaces and anchor in none of the structures (anchors = first 15% of the
% users in front of both, and of the users in front of each surface, in dataset order)
    n = height(P);  any_ = false(n,1);
    b = find(P.front1 & P.front2);  any_(b(1:round(0.15*numel(b)))) = true;
    for e = 1:2
        f = find(P.(sprintf('front%d',e)));  any_(f(1:round(0.15*numel(f)))) = true;
    end
    c = P.front1 & P.front2 & ~any_;
end

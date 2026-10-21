%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% POSITION ESTIMATION vs CODEBOOK SIZE (1x1 ... 5x5)
%
% Same experiment as Main_Multi_Slot_Fusion (position chart target, per-codebook
% semi-supervised t-SNE on Log-Euclidean covariance dissimilarities, then global /
% local fusion, 3200 users, 15% anchors, noiseless covariance + 1e-6 regularization),
% but instead of one 11x11 codebook it evaluates K x K codebooks for K = 1..5, for two
% phase-step ranges per EMS:
%   narrow : 0..40 deg/element (the range of the paper's 11x11 codebook)
%   wide   : K steps spread evenly over -180..180 deg/element (beams over all directions;
%            for lambda/2 spacing sin(theta_beam) = sin(theta_specular) + step/180)
% Each EMS configuration's chart is computed once and shared by every K that contains it.
% Fusion, best codebook: the repo's global_weighting / local_weighting /
% evaluate_best_codebook, unchanged.
%
% Speed-ups that give the same numbers as run_channel_charting (verified on the AoA
% branch): element-wise EMS phases instead of the 3600x3600 diagonal product, closed-form
% covariance (each realization is H times a scalar), logm once per user for the LED,
% O(N) rank look-ups in the TW/CT score.
%
% Environment: POS_WORKERS (default 24), POS_RANGES ("narrow,wide"), POS_KS ("1,2,3,4,5")
% Outputs: results_position_codebook/ (README, trend figures and CSVs per range)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();
nW     = str2double(envs('POS_WORKERS','24'));
ranges = strtrim(strsplit(envs('POS_RANGES','narrow,wide'), ','));
Ks     = str2double(strsplit(envs('POS_KS','1,2,3,4,5'), ','));
outRoot = fullfile(root, 'results_position_codebook');

wide = @(K) sort(mod((0:K-1)*360/K + 180, 360) - 180);
LEV.narrow = {0, [0 40], [0 20 40], [0 40/3 80/3 40], [0 10 20 30 40]};
LEV.wide   = arrayfun(wide, 1:5, 'UniformOutput', false);
TXT.narrow = 'phase steps 0..40 deg/element (the range of the paper''s 11x11 codebook)';
TXT.wide   = 'phase steps spread evenly over -180..180 deg/element (beams over all directions)';

%% ------------------------------------------------------------------------
% DATA (as Main_Multi_Slot_Fusion)
%% ------------------------------------------------------------------------
Data  = load_dataset(params);
state = prepare_data(Data, params);
clear Data
Len = state.Len;  nA = params.supervised_Len;
DistT = squareform(pdist(state.R_xy));
fprintf('Users %d, anchors %d\n', Len, nA);

%% ------------------------------------------------------------------------
% CHARTS: one per unique EMS configuration (phase steps in deg/element)
%% ------------------------------------------------------------------------
cfgs = zeros(0,2);
for r = ranges, for K = Ks
    [a,b] = ndgrid(LEV.(r{1}){K}, LEV.(r{1}){K});
    cfgs = [cfgs; a(:) b(:)]; %#ok<AGROW>
end, end
cfgs = unique(round(cfgs,6), 'rows', 'stable');
fprintf('%d unique EMS configurations to chart\n', size(cfgs,1));

pool = gcp('nocreate'); if isempty(pool), pool = parpool('Processes', nW); end
Hd = permute(double(state.H_dir), [2 3 1]);                 % Nb x Nr x U
X1 = reshape(permute(state.H_i_EMS1,[2 3 1]), size(state.H_i_EMS1,2), []);
X2 = reshape(permute(state.H_i_EMS2,[2 3 1]), size(state.H_i_EMS2,2), []);
Ho1 = double(state.H_o_EMS1);  Ho2 = double(state.H_o_EMS2);
Nb = size(Hd,1);  Nr = size(Hd,2);

CH = cell(size(cfgs,1),1);
futs = parallel.FevalFuture.empty;  fidx = [];  next = 1;  done = 0;
while next <= size(cfgs,1) || ~isempty(futs)
    while next <= size(cfgs,1) && numel(futs) < 2*nW
        D = led_input(cfgs(next,:), Hd, X1, X2, Ho1, Ho2, Nb, Nr, Len, params);
        futs(end+1) = parfeval(pool, @chart_job, 1, D, state.R_xy, nA, ...
                               params.tsne_perplexity, DistT, 0.05*Len); %#ok<SAGROW>
        fidx(end+1) = next;  next = next + 1; %#ok<SAGROW>
    end
    [i, out] = fetchNext(futs);
    CH{fidx(i)} = out;  futs(i) = [];  fidx(i) = [];  done = done + 1;
    if mod(done,10) == 0, fprintf('  %d/%d charts (%.1f min)\n', done, size(cfgs,1), toc(t0)/60); end
end
save(fullfile(root,'results_position_codebook_charts.mat'), 'cfgs', 'CH', '-v7.3');   % *.mat ignored

%% ------------------------------------------------------------------------
% EVALUATION PER RANGE AND K (fusion with the repo's functions)
%% ------------------------------------------------------------------------
unl = nA+1:Len;
ALL = {};
for r = ranges
    rDir = fullfile(outRoot, r{1});
    T = {};
    for K = Ks
        [a,b] = ndgrid(LEV.(r{1}){K}, LEV.(r{1}){K});
        [~, loc] = ismember(round([a(:) b(:)],6), cfgs, 'rows');
        C  = numel(loc);  R = [CH{loc}];
        rp = struct();
        rp.PosEstimates        = cat(3, R.xy);
        rp.ErrorPerUser        = [R.err];
        rp.ErrorPerCodebook    = arrayfun(@(o) mean(o.err(unl)), R).';
        rp.CT_PerUser          = [R.tc1];        % main's naming (column 1 = trustworthiness)
        rp.TW_PerUser          = [R.tc2];        % main's naming (column 2 = continuity)
        rp.CT_PerCodebook_mean = mean(rp.CT_PerUser(unl,:)).';
        rp.CT_PerCodebook_p90  = prctile(rp.CT_PerUser(unl,:),90).';
        rp.TW_PerCodebook_mean = mean(rp.TW_PerUser(unl,:)).';
        rp.TW_PerCodebook_p90  = prctile(rp.TW_PerUser(unl,:),90).';
        rg = global_weighting(rp, state, params);
        rl = local_weighting(rp, state, params);
        bb = evaluate_best_codebook(rp, state, params);

        M = {'global_fusion', rg.xy; 'local_fusion', rl.xy; 'best_single', bb.xy};
        rows = cell(0, 9);
        for m = 1:size(M,1)
            rows(end+1,:) = [M(m,1), num2cell(stats(M{m,2}, state.R_xy, DistT, unl))]; %#ok<AGROW>
        end
        S = zeros(C, 8);
        for c = 1:C, S(c,:) = stats(rp.PosEstimates(:,:,c), state.R_xy, DistT, unl); end
        rows(end+1,:) = [{'typical_single'}, num2cell(mean(S,1))];
        Tk = cell2table(rows, 'VariableNames', {'method','err_mean_m','err_median_m','err_p90_m', ...
                        'err_rmse_m','TW_mean','TW_p10','CT_mean','CT_p10'});
        Tk = [table(repmat(r(1),height(Tk),1), repmat(K,height(Tk),1), repmat(K^2,height(Tk),1), ...
              'VariableNames', {'range','K','acquisitions'}), Tk];
        kDir = fullfile(rDir, sprintf('K%d',K));  if ~exist(kDir,'dir'), mkdir(kDir); end
        writetable(Tk, fullfile(kDir,'summary.csv'));
        Tc = table((1:C).', round(a(:),2), round(b(:),2), rp.ErrorPerCodebook, ...
                   rp.TW_PerCodebook_mean, rp.CT_PerCodebook_mean, rg.weights(:), ...
            'VariableNames', {'config','step1_deg','step2_deg','err_mean_m', ...
                              'col2_continuity_mean','col1_trustworthiness_mean','global_weight'});
        writetable(Tc, fullfile(kDir,'per_configuration.csv'));
        plot_cdf(M, state.R_xy, unl, sprintf('%s codebook, %dx%d', r{1}, K, K), fullfile(kDir,'cdf.png'));
        T{end+1} = Tk; %#ok<SAGROW>
    end
    T = vertcat(T{:});
    writetable(T, fullfile(rDir,'trend.csv'));
    plot_trend(T, sprintf('%s codebook (%s)', r{1}, TXT.(r{1})), fullfile(rDir,'trend.png'));
    write_readme(rDir, r{1}, TXT.(r{1}), LEV.(r{1}), Ks, T, Len, nA);
    ALL{end+1} = T; %#ok<SAGROW>
end
ALL = vertcat(ALL{:});
writetable(ALL, fullfile(outRoot,'comparison.csv'));
plot_compare(ALL, fullfile(outRoot,'comparison.png'));
write_index(outRoot, ALL, ranges, TXT, toc(t0)/60, Len, nA);
fprintf('\n=== DONE (%.1f min) -> %s ===\n', toc(t0)/60, outRoot);


%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function D = led_input(cfg, Hd, X1, X2, Ho1, Ho2, Nb, Nr, Len, params)
% channel + covariance + Log-Euclidean dissimilarity of one EMS configuration
% (same computation as run_channel_charting + CompDismilarityLogEuclidean)
    th1 = exp(-1j*kron(deg2rad(cfg(1))*(1:params.M), ones(1,params.N))).';
    th2 = exp(-1j*kron(deg2rad(cfg(2))*(1:params.M), ones(1,params.N))).';
    H = Hd + reshape(Ho1*(th1.*X1), Nb, Nr, Len) + reshape(Ho2*(th2.*X2), Nb, Nr, Len);
    key = sprintf('%.6f_%.6f', cfg);
    rs  = RandStream('mt19937ar', 'Seed', mod(sum(double(key).*(1:numel(key))), 2^31-2) + 1);
    L = complex(zeros(Len, Nb*Nb));
    for u = 1:Len
        v  = randn(rs, 2, params.N_real);
        g2 = mean((v(1,:).^2 + v(2,:).^2)/2);
        Hu = H(:,:,u);
        R  = (Hu*Hu')*g2 + 1e-6*eye(Nb);
        Lu = logm(R);
        L(u,:) = Lu(:).';
    end
    D = squareform(pdist([real(L), imag(L)]));
end

function out = chart_job(D, xy, nA, perplexity, DistT, KK)
    out.xy  = tsne_sem(D, xy, 1:nA, perplexity);
    out.err = sqrt(sum((out.xy - xy).^2, 2));
    TC = fast_drscore(DistT, squareform(pdist(out.xy)), KK);
    out.tc1 = TC(:,1);  out.tc2 = TC(:,2);
end

function s = stats(xy, xyT, DistT, unl)
% position error [m] of the unlabeled users and their trustworthiness / continuity
% vs the true positions (K = 5% of the users, as compute_all_metrics)
    e  = sqrt(sum((xy(unl,:) - xyT(unl,:)).^2, 2));
    TC = fast_drscore(DistT, squareform(pdist(xy)), 0.05*size(xy,1));
    tw = TC(unl,1);  ct = TC(unl,2);
    s  = [mean(e), median(e), prctile(e,90), sqrt(mean(e.^2)), mean(tw), prctile(tw,10), mean(ct), prctile(ct,10)];
end

function TC = fast_drscore(D0, D, K)
% identical to drscore_pointwise(D0, D, K, 'TC') = [Trustworthiness Continuity]
    N = size(D0,1);  k = floor(K);
    TW = zeros(N,1); CT = zeros(N,1); rV = zeros(1,N); rU = zeros(1,N);
    for i = 1:N
        d0 = D0(i,:); d0(i) = 1e6;  d = D(i,:); d(i) = 1e6;
        [~,Vn] = sort(d0); [~,Un] = sort(d);
        rV(Vn) = 1:N; rU(Un) = 1:N;
        Vk = Vn(1:k); Uk = Un(1:k);
        inU = false(1,N); inU(Uk) = true;  inV = false(1,N); inV(Vk) = true;
        CT(i) = sum(rU(Vk(~inU(Vk))) - K);
        TW(i) = sum(rV(Uk(~inV(Uk))) - K);
    end
    c = 2/(K*(2*N-3*K-1));
    TC = [1 - c*TW, 1 - c*CT];
end

function plot_cdf(M, xyT, unl, ttl, file)
    f = figure('Position',[100 100 650 450]); hold on; grid on; box on;
    for m = 1:size(M,1)
        e = sqrt(sum((M{m,2}(unl,:) - xyT(unl,:)).^2, 2));
        plot(sort(e), linspace(0,1,numel(e)), 'LineWidth',1.6, 'DisplayName', strrep(M{m,1},'_',' '));
    end
    xlim([0 80]); xlabel('Position error [m]'); ylabel('CDF'); legend('Location','southeast'); title(ttl);
    exportgraphics(f, file, 'Resolution', 120); close(f);
end

function plot_trend(T, ttl, file)
    f = figure('Position',[100 100 1500 420]);
    tl = tiledlayout(f, 1, 4, 'TileSpacing','compact');
    Y = {'err_mean_m','mean position error [m]'; 'err_p90_m','90th percentile error [m]'; ...
         'TW_mean','trustworthiness (mean)'; 'CT_mean','continuity (mean)'};
    meths = unique(T.method, 'stable');
    for p = 1:4
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for m = meths.'
            r = sortrows(T(strcmp(T.method,m{1}),:), 'acquisitions');
            ls = '-o'; if contains(m{1},'single'), ls = ':s'; end
            plot(ax, r.acquisitions, r.(Y{p,1}), ls, 'LineWidth',1.5, 'DisplayName', strrep(m{1},'_',' '));
        end
        set(ax,'XScale','log'); xticks(ax, [1 4 9 16 25]);
        xlabel(ax,'EMS configurations (K^2)'); ylabel(ax, Y{p,2});
        if p == 1, legend(ax,'Location','northeast'); end
    end
    title(tl, ttl);
    exportgraphics(f, file, 'Resolution', 120); close(f);
end

function plot_compare(T, file)
    f = figure('Position',[100 100 1500 420]);
    tl = tiledlayout(f, 1, 4, 'TileSpacing','compact');
    Y = {'err_mean_m','mean position error [m]'; 'err_p90_m','90th percentile error [m]'; ...
         'TW_mean','trustworthiness (mean)'; 'CT_mean','continuity (mean)'};
    for p = 1:4
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for r = unique(T.range,'stable').'
            for m = {'local_fusion','typical_single'}
                q = sortrows(T(strcmp(T.range,r{1}) & strcmp(T.method,m{1}),:), 'acquisitions');
                ls = '-o'; if contains(m{1},'single'), ls = ':s'; end
                plot(ax, q.acquisitions, q.(Y{p,1}), ls, 'LineWidth',1.5, ...
                     'DisplayName', sprintf('%s, %s', r{1}, strrep(m{1},'_',' ')));
            end
        end
        set(ax,'XScale','log'); xticks(ax, [1 4 9 16 25]);
        xlabel(ax,'EMS configurations (K^2)'); ylabel(ax, Y{p,2});
        if p == 1, legend(ax,'Location','northeast'); end
    end
    title(tl, 'Position estimation vs codebook size: narrow vs wide phase-step range');
    exportgraphics(f, file, 'Resolution', 120); close(f);
end

function write_readme(rDir, name, txt, lev, Ks, T, Len, nA)
    L = {sprintf('# Position estimation vs codebook size: %s range', name), '', ...
         sprintf('Codebook: %s. Users %d, anchors %d (15%%), noiseless covariance with the 1e-6 regularization, as the paper experiment (Main_Multi_Slot_Fusion); only the codebook changes. Generated %s by `experiments/Main_Position_CodebookSize.m`.', ...
                 txt, Len, nA, char(datetime('now','Format','yyyy-MM-dd HH:mm'))), '', ...
         '| K | phase steps per EMS [deg/element] | configurations | local fusion mean / p90 error [m] | global fusion mean [m] | best single mean [m] | typical single mean [m] | local fusion TW / CT |', ...
         '|---|---|---|---|---|---|---|---|'};
    for K = Ks
        g = @(m,c) T{strcmp(T.method,m) & T.K==K, c};
        L{end+1} = sprintf('| %dx%d | %s | %d | %.1f / %.1f | %.1f | %.1f | %.1f | %.3f / %.3f |', K, K, ...
            mat2str(round(lev{K},2)), K^2, g('local_fusion','err_mean_m'), g('local_fusion','err_p90_m'), ...
            g('global_fusion','err_mean_m'), g('best_single','err_mean_m'), g('typical_single','err_mean_m'), ...
            g('local_fusion','TW_mean'), g('local_fusion','CT_mean')); %#ok<AGROW>
    end
    L = [L, {'', '## Definitions', ...
        '- Errors are over the unlabeled users. `typical_single` = average over the K^2 single-configuration charts (what one random configuration gives); `best_single` = the configuration with the lowest mean error (oracle, as `evaluate_best_codebook`); `global_fusion` / `local_fusion` = the repo''s `global_weighting` / `local_weighting` over the K^2 charts. Fusion gain = typical_single - fused error.', ...
        '- TW / CT = trustworthiness / continuity of the chart vs the true positions (K = 5% of the users), mean over the unlabeled users (p10 also in the CSVs). Note: the repo''s `compute_all_metrics` and `run_channel_charting` store these with swapped names (its `TW_*` is column 2 of `drscore_pointwise` = continuity); here they are named correctly. The fusion weights use the repo''s quality (column 2) unchanged.', ...
        '- Files: `trend.png` / `trend.csv` (error, TW, CT vs codebook size), `K<k>/summary.csv`, `K<k>/per_configuration.csv` (each single configuration), `K<k>/cdf.png`.'}];
    fid = fopen(fullfile(rDir,'README.md'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
end

function write_index(outRoot, T, ranges, TXT, runtime, Len, nA)
    L = {'# Position estimation vs codebook size (1x1 ... 5x5)', '', ...
         sprintf(['Paper position experiment (3200 users, %d anchors, per-configuration t-SNE + fusion, ' ...
                  'noiseless, 1e-6 regularization) with K x K codebooks instead of 11x11, for two ' ...
                  'phase-step ranges. Runtime %.0f min. See each range folder for details.'], nA, runtime), ''};
    for r = ranges, L{end+1} = sprintf('- [`%s`](%s/): %s', r{1}, r{1}, TXT.(r{1})); end %#ok<AGROW>
    L = [L, {'', '| range | K | local fusion mean error [m] | p90 [m] | typical single mean [m] | fusion gain [m] | local fusion TW | local fusion CT |', ...
             '|---|---|---|---|---|---|---|---|'}];
    for r = ranges
        for K = unique(T.K).'
            g = @(m,c) T{strcmp(T.range,r{1}) & strcmp(T.method,m) & T.K==K, c};
            L{end+1} = sprintf('| %s | %dx%d | %.1f | %.1f | %.1f | %.1f | %.3f | %.3f |', r{1}, K, K, ...
                g('local_fusion','err_mean_m'), g('local_fusion','err_p90_m'), g('typical_single','err_mean_m'), ...
                g('typical_single','err_mean_m') - g('local_fusion','err_mean_m'), ...
                g('local_fusion','TW_mean'), g('local_fusion','CT_mean')); %#ok<AGROW>
        end
    end
    L{end+1} = '';  L{end+1} = '`comparison.png`: narrow vs wide, local fusion and typical single configuration.';
    fid = fopen(fullfile(outRoot,'README.md'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
end

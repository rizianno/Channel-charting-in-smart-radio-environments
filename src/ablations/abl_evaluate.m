function T = abl_evaluate(scn, K, U, charts, R, CAT, noiseName, outDir)
% ABL_EVALUATE  Errors, CSVs and figures of one scenario at codebook size K x K.
%   R : containers.Map job key -> chart output (abl_chart)
% Returns the summary table (also written to outDir/summary.csv).
%
% Methods  F0: best_single (oracle), typical_single (mean over configurations), global_orig,
%              local_orig (fusion weights as the original code: quality vs true labels),
%              global_free, local_free (label-free quality vs input dissimilarities)
%          F1-F3: chart (single chart, no fusion)
% Subsets (per EMS e; anchors of the chart that estimates EMS e are never scored)
%   common  : users in front of both surfaces and not an anchor in any structure
%             (identical users for every scenario: the apples-to-apples comparison)
%   natural : all users this scenario estimates and that are in front of EMS e
%   reach   : natural and with a non-zero channel to EMS e
%   all     : (S3 only) every user, also behind EMS e, vs its geometric azimuth (worst case)
% Errors are estimate - true, wrapped to [-180, 180) deg; statistics use |error|.

    if ~exist(outDir,'dir'), mkdir(outDir); end
    jobs = abl_jobs(scn, K, charts, CAT, noiseName);
    isF0 = strcmp(scn.feature, 'F0');
    nCfg = numel(abl_levels(CAT, scn, K))^2;

    hasEst = false(U.n,2);  isAnch = false(U.n,2);
    EST = struct();
    if isF0
        mnames = {'global_orig','local_orig','global_free','local_free'};
        single = repmat({nan(U.n,2)}, 1, nCfg);
    else
        mnames = {'chart'};
    end
    for m = mnames, EST.(m{1}) = nan(U.n,2); end

    for c = charts
        u = c.users;  cols = c.cols;
        hasEst(u, cols) = true;
        isAnch(u(1:c.nA), cols) = true;
        jc = jobs(strcmp(arrayfun(@(j) j.chart.id, jobs, 'UniformOutput', false), c.id));
        if isF0
            outs = R(jc(1).key);
            for q = 2:numel(jc), outs(q) = R(jc(q).key); end
            E3 = cat(3, outs.est);
            [go, lo] = abl_fuse(E3, [outs.q_orig], c.nA, c.labels);
            [gf, lf] = abl_fuse(E3, [outs.q_free], c.nA, c.labels);
            EST.global_orig(u,cols) = go;  EST.local_orig(u,cols) = lo;
            EST.global_free(u,cols) = gf;  EST.local_free(u,cols) = lf;
            for k = 1:nCfg, single{k}(u,cols) = E3(:,:,k); end
        else
            o = R(jc.key);
            EST.chart(u,cols) = o.est;
        end
    end

    wrap = @(x) mod(x + 180, 360) - 180;
    jointChart = ~charts(1).perEMS;
    sub = struct();
    for e = 1:2
        nat = hasEst(:,e) & ~isAnch(:,e) & U.front(:,e);
        sub.common(:,e)  = U.both & ~U.anyAnchor & hasEst(:,e);
        sub.natural(:,e) = nat;
        sub.reach(:,e)   = nat & U.reach(:,e);
        if strcmp(scn.structure,'S3'), sub.all(:,e) = hasEst(:,e) & ~isAnch(:,e); end
    end
    subs = fieldnames(sub).';

    rows = {};
    addRows = @(name, Em) cellfun(@(s) stats_rows(name, s, wrap(Em - U.az), sub.(s), jointChart), ...
                                  subs, 'UniformOutput', false);
    for m = mnames
        r = addRows(m{1}, EST.(m{1}));  rows = [rows, r]; %#ok<AGROW>
    end
    T = vertcat(rows{:});
    if isF0
        Ts = cell(1,nCfg);
        for k = 1:nCfg
            r = addRows(sprintf('single_%d',k), single{k});  Ts{k} = vertcat(r{:});
        end
        T = [T; typical_and_best(Ts)];
    end
    T.scenario = repmat({scn.name}, height(T), 1);
    T.K = repmat(K, height(T), 1);
    T.acquisitions = repmat(abl_acquisitions(scn, K, jointChart), height(T), 1);
    T = movevars(T, {'scenario','K','acquisitions'}, 'Before', 1);
    writetable(T, fullfile(outDir,'summary.csv'));

    % per-user CSV of the main method
    if isF0, main = 'local_orig'; else, main = 'chart'; end
    Em = EST.(main);  Er = wrap(Em - U.az);
    P = table((1:U.n).', U.xy(:,1), U.xy(:,2), U.front(:,1), U.front(:,2), U.has_direct, ...
              U.reach(:,1), U.reach(:,2), isAnch(:,1), isAnch(:,2), U.az(:,1), U.az(:,2), ...
              Em(:,1), Em(:,2), Er(:,1), Er(:,2), ...
        'VariableNames', {'user','x_m','y_m','front1','front2','has_direct','reach1','reach2', ...
                          'anchor1','anchor2','az1_true_deg','az2_true_deg', ...
                          'az1_est_deg','az2_est_deg','err1_deg','err2_deg'});
    writetable(P, fullfile(outDir, sprintf('per_user_%s.csv', main)));

    % divergence check: a chart that blew up gives estimates far outside the physical range
    wf = fullfile(outDir, 'WARNING.txt');
    if isfile(wf), delete(wf); end
    out = abs(Em(hasEst)) > 180;
    if any(out)
        fid = fopen(wf, 'w');
        fprintf(fid, '%d of %d estimates (%.1f%%) are outside +-180 deg: the t-SNE chart diverged for them.\n', ...
                nnz(out), numel(out), 100*mean(out));
        fclose(fid);
    end

    plot_cdfs(EST, mnames, wrap, U, sub, jointChart, scn, K, fullfile(outDir,'cdf.png'));
    plot_maps(Er, sub.natural, U, 'abs', 30, sprintf('%s, K=%d, %s', scn.name, K, main), ...
              fullfile(outDir,'map_per_ems.png'));
    plot_maps(Er, sub.natural, U, 'signed', 30, sprintf('%s, K=%d, %s', scn.name, K, main), ...
              fullfile(outDir,'map_per_ems_signed.png'));
    if isfield(sub,'all')
        plot_maps(Er, sub.all, U, 'abs', 90, sprintf('%s, K=%d, %s, ALL users', scn.name, K, main), ...
                  fullfile(outDir,'map_per_ems_all_users.png'));
    end
end

% -------------------------------------------------------------------------
function T = stats_rows(method, subset, Err, S, jointChart)
    sets = {'EMS1', abs(Err(S(:,1),1)); 'EMS2', abs(Err(S(:,2),2)); ...
            'pooled', [abs(Err(S(:,1),1)); abs(Err(S(:,2),2))]};
    if jointChart
        b = S(:,1) & S(:,2);
        sets(end+1,:) = {'joint', sqrt(Err(b,1).^2 + Err(b,2).^2)};
    end
    n = size(sets,1);
    T = table(repmat({method},n,1), repmat({subset},n,1), sets(:,1), zeros(n,1), zeros(n,1), ...
              zeros(n,1), zeros(n,1), zeros(n,1), ...
              'VariableNames', {'method','subset','ems','n','median_deg','mean_deg','p90_deg','rmse_deg'});
    for i = 1:n
        x = sets{i,2};
        T{i,4:8} = [numel(x), median(x), mean(x), prctile(x,90), sqrt(mean(x.^2))];
    end
end

function T = typical_and_best(Ts)
% typical_single: mean over configurations of each statistic; best_single: the configuration
% with the lowest pooled mean error on the natural subset (oracle choice)
    T1 = Ts{1};
    v  = cellfun(@(t) t{:,4:8}, Ts, 'UniformOutput', false);
    T1{:,4:8} = mean(cat(3, v{:}), 3);
    T1.method(:) = {'typical_single'};
    crit = cellfun(@(t) t.mean_deg(strcmp(t.subset,'natural') & strcmp(t.ems,'pooled')), Ts);
    [~, kb] = min(crit);
    T2 = Ts{kb};  T2.method(:) = {sprintf('best_single(cfg%d)', kb)};
    T = [T1; T2];
end

function a = abl_acquisitions(scn, K, jointChart)
    a = K^2;
    if strcmp(scn.feature,'F3')
        if jointChart, a = a + 1; else, a = a + 2*K; end
    end
end

function plot_cdfs(EST, mnames, wrap, U, sub, jointChart, scn, K, file)
    np = 3 + jointChart;
    f = figure('Visible','off', 'Position',[100 100 420*np 380]);
    tl = tiledlayout(f, 1, np, 'TileSpacing','compact');
    ttl = {'EMS1','EMS2','pooled EMS1+EMS2','joint \surd(e_1^2+e_2^2)'};
    for p = 1:np
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for m = mnames
            Er = wrap(EST.(m{1}) - U.az);
            x  = pick(Er, sub.natural, p);
            plot(ax, sort(x), linspace(0,1,numel(x)), 'LineWidth',1.5, 'DisplayName', strrep(m{1},'_',' '));
            if isfield(sub,'all') && strcmp(m{1}, mnames{min(2,end)})
                xa = pick(Er, sub.all, p);
                plot(ax, sort(xa), linspace(0,1,numel(xa)), '--', 'LineWidth',1.2, ...
                     'DisplayName', [strrep(m{1},'_',' ') ', ALL users']);
            end
        end
        xlim(ax, [0 90]); xlabel(ax, '|AoA error| [deg]'); ylabel(ax, 'CDF'); title(ax, ttl{p});
        if p == 1, legend(ax, 'Location','southeast'); end
    end
    title(tl, strrep(sprintf('%s, K=%d (natural subset)', scn.name, K), '_', '\_'));
    exportgraphics(f, file, 'Resolution', 110); close(f);
end

function x = pick(Er, S, p)
    switch p
        case 1, x = abs(Er(S(:,1),1));
        case 2, x = abs(Er(S(:,2),2));
        case 3, x = [abs(Er(S(:,1),1)); abs(Er(S(:,2),2))];
        case 4, b = S(:,1) & S(:,2); x = sqrt(Er(b,1).^2 + Er(b,2).^2);
    end
end

function plot_maps(Er, S, U, mode, cmax, ttl, file)
    f = figure('Visible','off', 'Position',[100 100 1300 560]);
    tl = tiledlayout(f, 1, 2, 'TileSpacing','compact');
    for e = 1:2
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); axis(ax,'equal');
        s = S(:,e);
        if strcmp(mode,'signed')
            v = max(min(Er(s,e), cmax), -cmax); clim(ax, [-cmax cmax]); colormap(ax, bwr(256));
            lbl = sprintf('EMS%d error, estimate - true [deg]', e);
        else
            v = min(abs(Er(s,e)), cmax); clim(ax, [0 cmax]); colormap(ax, parula(256));
            lbl = sprintf('EMS%d |error| [deg] (>= %g saturates)', e, cmax);
        end
        scatter(ax, U.xy(s,1), U.xy(s,2), 10, v, 'filled');
        cb = colorbar(ax); ylabel(cb, lbl);
        o = 3 - e;
        plot(ax, U.EMS_position(o,1), U.EMS_position(o,2), 'ks', 'MarkerSize',8, 'MarkerFaceColor',[.8 .8 .8]);
        plot(ax, U.EMS_position(e,1), U.EMS_position(e,2), 'ks', 'MarkerSize',11, 'LineWidth',1.5, 'MarkerFaceColor','y');
        plot(ax, U.BS_position(1), U.BS_position(2), 'k^', 'MarkerSize',9, 'MarkerFaceColor','r');
        xlabel(ax,'x [m]'); ylabel(ax,'y [m]');
        title(ax, sprintf('Seen from EMS%d: median %.1f%s, n = %d', e, median(abs(Er(s,e))), char(176), nnz(s)));
    end
    title(tl, strrep(ttl,'_','\_'));
    exportgraphics(f, file, 'Resolution', 110); close(f);
end

function m = bwr(n)
    a = [0.23 0.30 0.75]; w = [0.97 0.97 0.97]; b = [0.71 0.02 0.15];
    h = floor(n/2); s = linspace(0,1,h)'; s2 = linspace(0,1,n-h)';
    m = [a + s.*(w-a); w + s2.*(b-w)];
end

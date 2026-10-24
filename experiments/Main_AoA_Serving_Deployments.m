%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMS DEPLOYMENTS COMPARED: EMS1 ALONE, EMS2 ALONE, BOTH (JOINT) - ONE EVALUATION FOR ALL SCHEMES
%
% Training deployments (what exists when the t-SNE training data are taken):
%   EMS1 alone / EMS2 alone : the other EMS does not exist (as Main_AoA_Serving_SingleEMS);
%   joint                   : both EMSs present, training slot k sets EMS1 to its k-th and EMS2
%                             to its k-th entry (as Main_AoA_Serving_NLoS_Size, informed/random/
%                             nested modes; charts reused from those caches).
% Evaluation views:
%   EMS1 / EMS2 : beam of that EMS only; SNR served by that EMS alone (+ the direct channel),
%                 users it reaches; angle = pointing error of the chosen beam (and, for t-SNE,
%                 the error of the estimated AoA itself);
%   joint       : (joint deployment only) both EMSs steered, SNR served by both, users reached by
%                 at least one EMS; angle errors pooled over both EMSs (every user-EMS pair).
% So "EMS1 alone" and "joint, EMS1 view" differ only in the training data.
% Schemes (sweeps up to SDP_SMAX; all incremental: each budget contains the previous one): greedy,
% hierarchical (full codebook and [clusters] = entries
% within 2 beamwidths of an anchor), exhaustive; t-SNE (per-user selection of the chart with the
% highest label-free trustworthiness) -> LUT, + R neighbours, + hierarchical warm-started at the
% LUT entry (full and [clusters]); training entries nested / informed / random (mean of draws).
% Populations (SDP_POP): fullnlos (charts from full-NLOS users), nolos (charts from every user
% without LoS; results split into full-NLOS and reflection-only users).
% Environment: SDP_POP ("fullnlos,nolos"), SDP_SIZES ("60r,60,120,180,240"), SDP_KS ("1:15"),
%   SDP_RS ("2,4,6"), SDP_SMAX (15), SDP_DRAWS (10), SDP_FOOT (0.15), SDP_NOISE (est_G60dB),
%   SDP_WORKERS (48), SDP_NEWCHARTS (1; 0 = only cached t-SNE charts, missing K skipped)
% Outputs: results_aoa/serving_deployments/<noise>/summary.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();  params.target = 'aoa';
CAT   = abl_catalog();
nz    = CAT.noise(strcmp({CAT.noise.name}, envs('SDP_NOISE','est_G60dB')));
pops  = strtrim(strsplit(envs('SDP_POP','fullnlos,nolos'), ','));
sizes = strtrim(strsplit(envs('SDP_SIZES','60r,60,120,180,240'), ','));
Ks    = eval(envs('SDP_KS','1:15'));  Kmax = max(Ks);
Rs    = str2double(strsplit(envs('SDP_RS','2,4,6'), ','));
Smax  = str2double(envs('SDP_SMAX','15'));
nDraw = str2double(envs('SDP_DRAWS','10'));
foot  = str2double(envs('SDP_FOOT','0.15'));
nW    = str2double(envs('SDP_WORKERS','48'));
newCh = str2double(envs('SDP_NEWCHARTS','1'));                  % 0: only cached t-SNE charts (missing K skipped)
k0 =2*pi/(3e8/30e9);  d = (3e8/30e9)/2;  lamc = 3e8/30e9;

Data = load_dataset(params);
U    = abl_users(Data, params);
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);  M0 = params.M;
Pd = squeeze(sum(sum(abs(Data.H_LOS_1bounce).^2,2),3)) / (Nb*Nr);
dd = sqrt(sum(([U.xy, params.UE_height*ones(U.n,1)] - U.BS_position).^2, 2));
los = U.has_direct & 10*log10(Pd ./ (lamc./(4*pi*dd)).^2) > -36;
rho = 10^((nz.tx_dBm - nz.noise_dBm)/10);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);  specular = -azBS;

pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', nW); end
outDir = fullfile(root,'results_aoa','serving_deployments', nz.name);
if ~exist(outDir,'dir'), mkdir(outDir); end
rows = {};
for pp = pops
    pop = pp{1};
    if strcmp(pop,'nolos'), popU = find(~los); else, popU = find(~U.has_direct); end   % order as the earlier runs
    nU = numel(popU);  isNone = ~U.has_direct(popU);
    reach = [any(Data.H_I_RIS1(popU,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(popU,:,:) ~= 0, [2 3])] & U.front(popU,:);
    Hd  = reshape(permute(double(Data.H_LOS_1bounce(popU,:,:)),[2 3 1]), Nb*Nr, nU);
    xyz = [U.xy(popU,:), params.UE_height*ones(nU,1)];
    ch = struct('users',{}, 'nA',{}, 'lab',{}, 'perp',{}, 'uA',{});
    isAnch = false(nU,2);
    for e = 1:2
        m = find(reach(:,e));  nA = round(params.Supervision_percentage * numel(m));
        ch(e).users = m;  ch(e).nA = nA;  ch(e).lab = U.az(popU(m), e);
        ch(e).perp = min(params.tsne_perplexity, round(numel(m)/5));  ch(e).uA = sind(ch(e).lab(1:nA));
        isAnch(m(1:nA), e) = true;
    end
    fit = fit_model(Data, U, popU, xyz, M0, d, k0);
    fprintf('\n=== population %s: %d users; EMS1 chart %d (%d anchors), EMS2 chart %d (%d)\n', pop, nU, ...
            numel(ch(1).users), ch(1).nA, numel(ch(2).users), ch(2).nA);

    for sz = sizes
        real_ = endsWith(sz{1}, 'r');  M = str2double(erase(sz{1}, 'r'));
        tag = sprintf('M%d%s', M, repmat('_raytraced', 1, real_));
        ctag = tag;  if strcmp(pop,'nolos'), ctag = ['nolos_' tag]; end
        wu = 0.886 * 2 / M;  N = ceil(2/wu);  N = N + mod(N,2);
        step = ((1:N) - N/2) * 360/N;
        A = ems_paths(Data, U, popU, xyz, fit, real_, M, M0, N, step, d, k0, Nb, Nr);   % {e}: (Nb*Nr) x nU x N
        G = cellfun(@(Ae) gains(Ae, Nb, Nr), A, 'UniformOutput', false);               % own-path gains nU x N
        lmax = @(H) lam_max(H, Nb, Nr);
        ub = arrayfun(@(e) mod(sind(specular(e)) + step/180 + 1, 2) - 1, 1:2, 'UniformOutput', false);
        P = struct('Sup',{}, 'sq',{}, 'lists',{});
        for e = 1:2
            near = abs(ch(e).uA - ub{e});
            P(e).Sup = find(any(near <= 2*wu, 1));
            [~, o] = sort(-sum(near(:,P(e).Sup) <= wu, 1));  P(e).sq = [P(e).Sup(o), setdiff(1:N, P(e).Sup)];
            no = nest_order(N);
            P(e).lists = {'nested', {no(1:Kmax)}; ...
                          'informed', {informed_order(ch(e).uA, ub{e}, foot, Kmax, max(2, round(0.02*ch(e).nA)))}; ...
                          'random', arrayfun(@(dI) randperm(RandStream('mt19937ar','Seed',1000*dI + e), N, Kmax), 1:nDraw, 'UniformOutput', false)};
        end
        fprintf('[%s] N = %d\n', tag, N);

        % ---- views: users, base SNR, groups
        scoredE = arrayfun(@(e) reach(:,e) & ~isAnch(:,e), 1:2, 'UniformOutput', false);
        base = {rho * lmax(Hd + pick(A{1}, N/2*ones(nU,1))), rho * lmax(Hd + pick(A{2}, N/2*ones(nU,1)))};
        baseJ = rho * lmax(Hd + pick(A{1}, N/2*ones(nU,1)) + pick(A{2}, N/2*ones(nU,1)));
        snrE = @(e, i) rho * lmax(Hd + pick(A{e}, i));
        snrJ = @(i1, i2) rho * lmax(Hd + pick(A{1}, i1) + pick(A{2}, i2));
        azTrue = U.az(popU, :);
        beamAz = @(e, i) asind(ub{e}(i(:))).';
        grpOf = @(m) {'full NLOS', m & isNone, isNone; 'reflection only', m & ~isNone, ~isNone};

        % ---- charts: EMS alone (per entry) and joint (per pair), from caches or computed
        Ra = cell(1,2);
        for e = 1:2
            ent = unique(cell2mat(cellfun(@(c) cell2mat(c), P(e).lists(:,2).', 'UniformOutput', false)));
            cf = fullfile(root,'cache_ablations', sprintf('single_%s_e%d_%s_%s.mat', nz.name, e, pop, tag));
            if ~newCh, Ra{e} = load_cache(cf); continue; end
            Ra{e} = alone_charts(load_cache(cf), cf, ent, A{e}(:, ch(e).users, :), Hd(:, ch(e).users), ch(e), ...
                                 sprintf('%s|single|e%d|%s', tag, e, pop), nz, params, pool, nW, Nb, Nr);
        end
        cfJ = fullfile(root,'cache_ablations', sprintf('pairs_%s_%s.mat', nz.name, ctag));
        Rj = load_cache(cfJ);
        for f = {'nested','nestedgrid'}
            Rn = load_cache(fullfile(root,'cache_ablations', sprintf('%s_%s_%s.mat', f{1}, nz.name, ctag)));
            for k = Rn.keys, if ~isKey(Rj, k{1}), Rj(k{1}) = Rn(k{1}); end, end
        end
        Pall = zeros(0,2);
        for t = 1:3
            for dI = 1:numel(P(1).lists{t,2}), Pall = [Pall; P(1).lists{t,2}{dI}(:), P(2).lists{t,2}{dI}(:)]; end %#ok<AGROW>
        end
        if newCh, Rj = joint_charts(Rj, cfJ, unique(Pall,'rows','stable'), A, Hd, ch, nz, ctag, params, pool, nW, Nb, Nr, nU); end
        fprintf('  charts ready (%.1f min)\n', toc(t0)/60);

        % ---- searches (independent of the training deployment): per EMS, then the joint view
        hdr = {pop, tag, M, N};
        add = @(rows, dep, view, train, scheme, K, R, s, snr, b, m, angB, angT) ...
              [rows; collect2([hdr, {dep, view, train, scheme, K, R, s}], snr, b, grpOf(m), angB, angT)];
        exh = cell(1,2);  exc = cell(1,2);
        for e = 1:2
            exh{e} = best1(G{e}, repmat(1:N, nU, 1));  exc{e} = best1(G{e}, repmat(P(e).Sup, nU, 1));
        end
        srch = {'exhaustive', @(e, s) exh{e}; 'exhaustive[clusters]', @(e, s) exc{e}; ...
                'greedy', @(e, s) best1(G{e}, repmat(1:min(s,N), nU, 1)); 'greedy[clusters]', @(e, s) best1(G{e}, repmat(P(e).sq(1:min(s,N)), nU, 1)); ...
                'hierarchical', @(e, s) hier1(G{e}, N, s); 'hierarchical[clusters]', @(e, s) P(e).Sup(hier1(G{e}(:,P(e).Sup), numel(P(e).Sup), s)).'};
        okJ = any(reach,2) & ~any(isAnch,2);
        for q = 1:size(srch,1)
            sw = 1:Smax;  if startsWith(srch{q,1},'exhaustive'), sw = 0; end
            for s = sw
                ii = {srch{q,2}(1, s), srch{q,2}(2, s)};  aeE = cell(1,2);  mE = cell(1,2);
                for e = 1:2
                    m = scoredE{e} & snrE(e, exh{e}).' > 0;
                    ae = abs(mod(beamAz(e, ii{e}) - azTrue(:,e) + 180, 360) - 180);
                    aeE{e} = ae;  mE{e} = m;
                    for dep = {sprintf('EMS%d alone', e), 'joint'}
                        rows = add(rows, dep{1}, sprintf('EMS%d', e), '-', srch{q,1}, 0, 0, s, snrE(e, ii{e}), base{e}, m, ae, nan(nU,1));
                    end
                end
                rows = add(rows, 'joint', 'joint', '-', srch{q,1}, 0, 0, s, snrJ(ii{1}, ii{2}), baseJ, okJ, ...
                           {aeE{1}, mE{1}, aeE{2}, mE{2}}, nan(nU,1));                   % angles pooled over both EMSs
            end
        end
        lutT = [lut1(azTrue(:,1), specular(1), step), lut1(azTrue(:,2), specular(2), step)];
        for e = 1:2
            m = scoredE{e} & snrE(e, exh{e}).' > 0;
            for dep = {sprintf('EMS%d alone', e), 'joint'}
                rows = add(rows, dep{1}, sprintf('EMS%d', e), '-', 'true AoA + LUT', 0, 0, 0, snrE(e, lutT(:,e)), base{e}, m, ...
                           abs(mod(beamAz(e, lutT(:,e)) - azTrue(:,e) + 180, 360) - 180), nan(nU,1));
            end
        end

        % ---- t-SNE: per training deployment, training list, K (random: mean over draws)
        for dep = {'EMS1 alone', 'EMS2 alone', 'joint'}
            for t = 1:3
                tr = P(1).lists{t,1};
                for K = Ks
                    rr = {};
                    for dI = 1:numel(P(1).lists{t,2})
                        az = nan(nU,2);  ee = 1:2;
                        if ~strcmp(dep{1}, 'joint'), ee = sscanf(dep{1}, 'EMS%d'); end
                        have = true;                                % all K charts cached?
                        for e = ee
                            li = P(e).lists{t,2}{dI}(1:K);
                            if strcmp(dep{1}, 'joint')
                                pr = zeros(K,2);  pr(:,e) = li;  pr(:,3-e) = P(3-e).lists{t,2}{dI}(1:K);
                                have = have && all(arrayfun(@(q) isKey(Rj, nkey(nz.name, ctag, e, pr(q,:))), 1:K));
                            else
                                have = have && all(arrayfun(@(i) isKey(Ra{e}, sprintf('%d', i)), li));
                            end
                        end
                        if ~have, continue; end
                        for e = ee
                            li = P(e).lists{t,2}{dI}(1:K);
                            if strcmp(dep{1}, 'joint')
                                lo = P(3-e).lists{t,2}{dI}(1:K);  pr = zeros(K,2);  pr(:,e) = li;  pr(:,3-e) = lo;
                                outs = Rj(nkey(nz.name, ctag, e, pr(1,:)));
                                for q = 2:K, outs(q) = Rj(nkey(nz.name, ctag, e, pr(q,:))); end
                            else
                                outs = Ra{e}(sprintf('%d', li(1)));
                                for q = 2:K, outs(q) = Ra{e}(sprintf('%d', li(q))); end
                            end
                            az(ch(e).users, e) = select_chart([outs.est], [outs.tw], ch(e).nA, ch(e).lab);
                        end
                        lut = [lut1(az(:,1), specular(1), step), lut1(az(:,2), specular(2), step)];
                        lpos = cell(1,2);
                        for e = 1:2, [~, lpos{e}] = min(abs(ub{e}(lut(:,e)).' - ub{e}(P(e).Sup)), [], 2); end
                        sch = {'tsne', 0, K, @(e) lut(:,e)};
                        for R = Rs
                            sch(end+1,:) = {'tsne+refine', R, K+R+1, @(e) best1(G{e}, nbr(lut(:,e), R, N))}; %#ok<AGROW>
                            sch(end+1,:) = {'tsne+refine[clusters]', R, K+R+1, @(e) best1(G{e}, P(e).Sup(nbr(lpos{e}, R, numel(P(e).Sup))))}; %#ok<AGROW>
                        end
                        for b = 1:(Smax - K)
                            sch(end+1,:) = {'tsne+hierarchical', 0, K+b, @(e) warm1(G{e}, lut(:,e), N, b, round(N/8))}; %#ok<AGROW>
                            sch(end+1,:) = {'tsne+hierarchical[clusters]', 0, K+b, @(e) P(e).Sup(warm1(G{e}(:,P(e).Sup), lpos{e}, numel(P(e).Sup), b, max(1, round(numel(P(e).Sup)/8)))).'}; %#ok<AGROW>
                        end
                        for q = 1:size(sch,1)
                            if sch{q,3} > Smax, continue; end
                            ii = {sch{q,4}(1), sch{q,4}(2)};  aBe = cell(1,2);  aTe = cell(1,2);  mE = cell(1,2);
                            for e = ee
                                m = scoredE{e} & snrE(e, exh{e}).' > 0;
                                aB = abs(mod(beamAz(e, ii{e}) - azTrue(:,e) + 180, 360) - 180);
                                aT = abs(mod(az(:,e) - azTrue(:,e) + 180, 360) - 180);
                                aBe{e} = aB;  aTe{e} = aT;  mE{e} = m;
                                rr = add(rr, dep{1}, sprintf('EMS%d', e), tr, sch{q,1}, K, sch{q,2}, sch{q,3}, snrE(e, ii{e}), base{e}, m, aB, aT);
                            end
                            if strcmp(dep{1}, 'joint')                     % angles pooled over both EMSs
                                rr = add(rr, 'joint', 'joint', tr, sch{q,1}, K, sch{q,2}, sch{q,3}, snrJ(ii{1}, ii{2}), baseJ, okJ, ...
                                         {aBe{1}, mE{1}, aBe{2}, mE{2}}, {aTe{1}, mE{1}, aTe{2}, mE{2}});
                            end
                        end
                    end
                    rows = [rows; avg_rows(rr, 14:25, 5:11)]; %#ok<AGROW>
                end
            end
            fprintf('  %s done (%.1f min)\n', dep{1}, toc(t0)/60);
        end
        clear A G
        write_table(outDir, rows);
    end
end
fprintf('\n=== DONE (%.1f min) ===\n', toc(t0)/60);

%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function write_table(outDir, rows)
    T = cell2table(rows, 'VariableNames', {'population','surface','M','N','deployment','view','training','scheme','K','R','sweeps', ...
        'group','n_users','dsnr_median_dB','dsnr_mean_dB','dsnr_p10_dB','snr_median_dB','beam_err_median_deg','beam_err_rmse_deg', ...
        'beam_err_p90_deg','tsne_err_median_deg','tsne_err_rmse_deg','tsne_err_p90_deg','se_mean','dse_mean','n_draws','x'});
    T.x = [];
    writetable(T, fullfile(outDir,'summary.csv'));
end

function r = collect2(key, snr, base, grp, angB, angT)
% one row per group: key (11 fields), group, n, SNR gain median/mean/p10, SNR median, beam and
% t-SNE angle error median/RMSE/p90, mean SE and its gain, 1 (draw count), placeholder.
% angB / angT: per-user errors (masked by the group), or {err1, mask1, err2, mask2} = pooled over
% both EMSs (every user-EMS pair of the group, each with its own EMS mask)
    r = {};
    pool2 = @(a, G) [a{1}(a{2} & G); a{3}(a{4} & G)];
    for gi = 1:size(grp,1)
        m = grp{gi,2};  if ~any(m), continue; end
        x = snr(m);  b = base(m);  dd = 10*log10(x(:) ./ b(:));  dd = dd(isfinite(dd));
        if iscell(angB), aB = pool2(angB, grp{gi,3}); else, aB = angB(m); end
        if iscell(angT), aT = pool2(angT, grp{gi,3}); else, aT = angT(m); end
        r(end+1,:) = [key, {grp{gi,1}, nnz(m), median(dd), mean(dd), prctile(dd,10), median(10*log10(x)), ...
            median(aB), sqrt(mean(aB.^2)), prctile(aB,90), median(aT), sqrt(mean(aT.^2)), prctile(aT,90), ...
            mean(log2(1+x)), mean(log2(1+x) - log2(1+b)), 1, 0}]; %#ok<AGROW>
    end
end

function out = avg_rows(rr, numCols, keyCols)
% rows of several random draws -> one row per key (key columns + group), numeric columns averaged
    if isempty(rr), out = rr; return; end
    key = arrayfun(@(i) strjoin(cellfun(@num2str, rr(i,[keyCols 12]), 'UniformOutput', false), '|'), (1:size(rr,1)).', 'UniformOutput', false);
    [u, ia, ic] = unique(key, 'stable');
    if numel(u) == size(rr,1), out = rr; return; end
    out = rr(ia,:);  X = cell2mat(rr(:,numCols));
    for j = 1:numel(u), out(j,numCols) = num2cell(mean(X(ic==j,:), 1, 'omitnan')); out{j,26} = nnz(ic==j); end
end

function fit = fit_model(Data, U, popU, xyz, M0, d, k0)
% geometric model per EMS fitted on the 60 x 60 data (as Main_AoA_Serving_NLoS_Size)
    fit = struct();
    for e = 1:2
        Ho = double(Data.(sprintf('H_o_RIS%d',e)));  [Uo, ~, Vo] = svds(Ho, 1);
        Hi = Data.(sprintf('H_I_RIS%d',e));  Hi = Hi(popU,:,:);
        best = -inf;
        for sm = [1 -1], for sn = [1 -1], for pc = [1 -1]
            b0 = sph(elem_pos(U.EMS_position(e,:), M0, d, sm, sn), U.BS_position, k0, pc);
            cb = abs(Vo' * conj(b0)) / norm(b0);
            if cb > best, best = cb; ori = [sm sn pc]; end
        end, end, end
        P0 = elem_pos(U.EMS_position(e,:), M0, d, ori(1), ori(2));
        b0 = sph(P0, U.BS_position, k0, ori(3));
        S0 = sph(P0, xyz, k0, ori(3));  alpha = zeros(numel(popU), size(Hi,3));
        for u = 1:numel(popU), s0 = S0(:,u); alpha(u,:) = (s0' * squeeze(Hi(u,:,:))) / (s0' * s0); end
        fit(e).ori = ori;  fit(e).a = Uo;  fit(e).gamma = (Uo' * Ho * conj(b0)) / (b0.' * conj(b0));  fit(e).alpha = alpha;
    end
end

function A = ems_paths(Data, U, popU, xyz, fit, real_, M, M0, N, step, d, k0, Nb, Nr)
% per EMS, per entry: the channel through that EMS of every user, (Nb*Nr) x nU x N
    nU = numel(popU);  A = cell(1,2);
    for e = 1:2
        Ce = zeros(N, Nb*Nr, nU);
        if real_
            Hi = Data.(sprintf('H_I_RIS%d',e));  Hi = Hi(popU,:,:);
            Ho = double(Data.(sprintf('H_o_RIS%d',e)));  X = reshape(permute(Hi,[2 3 1]), M0^2, []);
            for i = 1:N
                th = exp(-1j*kron(deg2rad(step(i))*(1:M0), ones(1,M0))).';
                Ce(i,:,:) = reshape(reshape(Ho*(th.*X), Nb, Nr, nU), Nb*Nr, nU);
            end
        else
            P = elem_pos(U.EMS_position(e,:), M, d, fit(e).ori(1), fit(e).ori(2));
            bB = sph(P, U.BS_position, k0, fit(e).ori(3));  Y = zeros(M, nU);
            for mm = 1:M, k = (mm-1)*M + (1:M); Y(mm,:) = bB(k).' * sph(P(k,:), xyz, k0, fit(e).ori(3)); end
            c = exp(-1j * deg2rad(step(:)) * (1:M)) * Y;
            for r = 1:Nr
                Ce(:, (r-1)*Nb + (1:Nb), :) = reshape(fit(e).gamma*fit(e).a, 1, Nb, 1) .* reshape(c .* fit(e).alpha(:,r).', N, 1, nU);
            end
        end
        A{e} = permute(Ce, [2 3 1]);
    end
end

function Rmap = alone_charts(Rmap, cacheFile, ent, Ae, Hde, che, skey, nz, params, pool, nW, Nb, Nr)
% one chart per entry with only this EMS deployed (as Main_AoA_Serving_SingleEMS; same seeds)
    need = ent(~arrayfun(@(i) isKey(Rmap, sprintf('%d', i)), ent));
    futs = parallel.FevalFuture.empty;  fk = {};
    for i = need
        F = logcov_features(reshape(Hde + Ae(:,:,i), Nb, Nr, []), sprintf('%s|%d', skey, i), nz, params);
        futs(end+1) = parfeval(pool, @srv_chart_twct, 1, squareform(pdist(F)), che.lab, che.nA, che.perp, params.K_fraction); %#ok<AGROW>
        fk{end+1} = sprintf('%d', i); %#ok<AGROW>
        while numel(futs) >= 2*nW, [j, out] = fetchNext(futs); Rmap(fk{j}) = out; futs(j) = []; fk(j) = []; end
    end
    while ~isempty(futs), [j, out] = fetchNext(futs); Rmap(fk{j}) = out; futs(j) = []; fk(j) = []; end
    if ~isempty(need), save_cache(cacheFile, Rmap); end
    fprintf('  alone charts: %d computed\n', numel(need));
end

function Rmap = joint_charts(Rmap, cacheFile, P, A, Hd, ch, nz, ctag, params, pool, nW, Nb, Nr, nU)
% one chart per (pair, EMS) with both EMSs present (as Main_AoA_Serving_NLoS_Size; same seeds/keys)
    futs = parallel.FevalFuture.empty;  fk = {};  nCh = 0;
    for q = 1:size(P,1)
        pr = P(q,:);
        need = arrayfun(@(e) ~isKey(Rmap, nkey(nz.name, ctag, e, pr)), 1:2);
        if ~any(need), continue; end
        H = reshape(Hd + A{1}(:,:,pr(1)) + A{2}(:,:,pr(2)), Nb, Nr, nU);
        if pr(1) == pr(2), skey = sprintf('%s|%d', ctag, pr(1)); else, skey = sprintf('%s|%d_%d', ctag, pr); end
        F = logcov_features(H, skey, nz, params);
        for e = find(need)
            futs(end+1) = parfeval(pool, @srv_chart_twct, 1, squareform(pdist(F(ch(e).users,:))), ch(e).lab, ch(e).nA, ch(e).perp, params.K_fraction); %#ok<AGROW>
            fk{end+1} = nkey(nz.name, ctag, e, pr); %#ok<AGROW>
        end
        while numel(futs) >= 2*nW, [j, out] = fetchNext(futs); Rmap(fk{j}) = out; futs(j) = []; fk(j) = []; nCh = nCh + 1; end
    end
    while ~isempty(futs)
        [j, out] = fetchNext(futs);  Rmap(fk{j}) = out;  futs(j) = [];  fk(j) = [];  nCh = nCh + 1;
        if mod(nCh, 100) == 0, save_cache(cacheFile, Rmap); end
    end
    if nCh > 0, save_cache(cacheFile, Rmap); end
    fprintf('  joint charts: %d computed\n', nCh);
end

function k = nkey(nzName, ctag, e, pr)
    if pr(1) == pr(2), k = sprintf('%s|%s|nest|e%d|%d', nzName, ctag, e, pr(1));
    else, k = sprintf('%s|%s|ngrid|e%d|%d_%d', nzName, ctag, e, pr(1), pr(2)); end
end

function x = select_chart(E, QT, nA, lab)
% per user, the chart with the highest label-free trustworthiness (anchors keep their labels)
    [n, K] = size(E);  [~, k] = max(QT, [], 2);  x = E(sub2ind([n K], (1:n).', k));  x(1:nA) = lab(1:nA);
end

function P = elem_pos(c, M, d, sm, sn)
    o = ((1:M) - (M+1)/2) * d;
    [n_, m_] = meshgrid(o, o);  m_ = m_.';  n_ = n_.';
    P = c + [zeros(M^2,1), sm*m_(:), sn*n_(:)];
end

function S = sph(P, Q, k0, pc)
    r = sqrt((P(:,1) - Q(:,1).').^2 + (P(:,2) - Q(:,2).').^2 + (P(:,3) - Q(:,3).').^2);
    S = exp(-1j*pc*k0*r) ./ r;
end

function F = logcov_features(H, key, noise, params)
% log-covariance features (same model and seeds as the earlier runs)
    [Nb, Nr, nU] = size(H);
    key = [key '|' noise.name];
    rs  = RandStream('mt19937ar', 'Seed', mod(sum(double(key).*(1:numel(key))), 2^31-2) + 1);
    sn2 = 10^((noise.noise_dBm - noise.tx_dBm)/10);
    noisy = ~isnan(noise.pilot_gain_dB);
    if noisy, se2 = sn2 / 10^(noise.pilot_gain_dB/10); end
    flo = 2*sn2;  if isfield(noise,'floor_abs') && ~isnan(noise.floor_abs), flo = noise.floor_abs; end
    Nre = params.N_real;  iu = find(triu(true(Nb),1));
    F = zeros(nU, Nb + 2*numel(iu));
    for u = 1:nU
        Hu = H(:,:,u);  v = randn(rs, 2*Nre, 1);  g2 = sum(v.^2)/2;
        if noisy
            Am = sqrt(g2)*Hu + sqrt(se2/2)*complex(randn(rs,Nb,Nr), randn(rs,Nb,Nr));
            Rm = (Am*Am' + se2*cwishart(rs, Nb, Nr*(Nre-1))) / Nre;
        else
            Rm = (g2/Nre)*(Hu*Hu') + flo*eye(Nb);
        end
        Rm = (Rm + Rm')/2;  [V, ev] = eig(Rm, 'vector');  L = (V .* log(ev).') * V';
        F(u,:) = [real(diag(L)).', sqrt(2)*real(L(iu)).', sqrt(2)*imag(L(iu)).'];
    end
end

function W = cwishart(rs, p, n)
    T = zeros(p);
    for i = 1:p, T(i,i) = sqrt(sum(randn(rs, 2*(n-i+1), 1).^2)/2); end
    lo = find(tril(true(p),-1));  T = complex(T);
    T(lo) = complex(randn(rs,numel(lo),1), randn(rs,numel(lo),1))/sqrt(2);
    W = T*T';
end

function l = lam_max(H, Nb, Nr)
    H  = reshape(H, Nb, Nr, []);
    h1 = squeeze(H(:,1,:));  h2 = squeeze(H(:,2,:));
    if isvector(h1), h1 = h1(:); h2 = h2(:); end
    g11 = sum(abs(h1).^2,1);  g22 = sum(abs(h2).^2,1);  g12 = sum(conj(h1).*h2,1);
    l = (g11+g22)/2 + sqrt(((g11-g22)/2).^2 + abs(g12).^2);
end

function G = gains(Ae, Nb, Nr)
    G = zeros(size(Ae,2), size(Ae,3));
    for i = 1:size(Ae,3), G(:,i) = lam_max(Ae(:,:,i), Nb, Nr).'; end
end

function H = pick(Ae, idx)
    [dd, n, ~] = size(Ae);
    H = reshape(Ae(:, (idx(:).'-1)*n + (1:n)), dd, n);
end

function b = best1(g, c)
    n = size(c,1);  [~, k] = max(g((c - 1)*n + (1:n).'), [], 2);  b = c((k - 1)*n + (1:n).');
end

function c = nbr(b, R, N)
    c = mod(b(:) - 1 + [0, -(R/2):-1, 1:(R/2)], N) + 1;
end

function b = hier1(g, N, s)
% hierarchical search, incremental (design B): the first min(s,4) probes are the first entries of
% one fixed ever-finer order (N/2, N, N/4, 3N/4: uniform on the circle once all four are taken),
% then +-delta around the best so far, delta halved; every budget contains the previous one
    n = size(g,1);  o = nest_order(N);
    b = best1(g, repmat(o(1:min(s,4)), n, 1));  used = min(s,4);  delta = N/4;
    while used < s && delta > 1
        delta = ceil(delta/2);
        for sgn = [-1 1]
            if used >= s, break; end
            cand = mod(b - 1 + sgn*delta, N) + 1;
            better = g((cand-1)*n + (1:n).') > g((b-1)*n + (1:n).');
            b(better) = cand(better);  used = used + 1;
        end
    end
end

function b = warm1(g, b, N, budget, delta0)
    n = size(g,1);  used = 1;  delta = delta0;
    while used < budget && delta >= 1
        for sgn = [-1 1]
            if used >= budget, break; end
            cand = mod(b - 1 + sgn*delta, N) + 1;
            better = g((cand-1)*n + (1:n).') > g((b-1)*n + (1:n).');
            b(better) = cand(better);  used = used + 1;
        end
        if delta == 1, break; end
        delta = ceil(delta/2);
    end
end

function sel = lut1(az, spec, step)
    want = mod(180*(sind(az) - sind(spec)) + 180, 360) - 180;
    [~, sel] = min(abs(mod(want - step + 180, 360) - 180), [], 2);
    sel(isnan(az)) = numel(step)/2;
end

function o = nest_order(N)
    f = [1/2 1];
    for lev = 2:12
        f = [f, (1:2:2^lev-1)/2^lev]; %#ok<AGROW>
        if numel(unique(round(f*N))) >= N, break; end
    end
    o = unique(min(max(round(f*N),1),N), 'stable');
end

function o = informed_order(uA, ub, F, Kmax, minA)
    o = [];  r = F;
    while numel(o) < Kmax && r > 1e-3
        cov = false(numel(uA),1);
        if ~isempty(o), cov = any(abs(uA - ub(o)) <= r, 2); end
        while numel(o) < Kmax
            unc = find(~cov);
            if numel(unc) < minA, break; end
            dd = abs(uA(unc) - ub);  inR = dd <= r;  cnt = sum(inR,1);  cnt(o) = -1;
            if max(cnt) < minA, break; end
            cand = find(cnt == max(cnt));
            md = arrayfun(@(c) mean(dd(inR(:,c),c)), cand);
            [~, j] = min(md);  c = cand(j);
            o(end+1) = c;  cov(unc(inR(:,c))) = true; %#ok<AGROW>
        end
        r = r/2;
    end
    if numel(o) < Kmax
        [~, k] = sort(min(abs(uA - ub), [], 1));  k = setdiff(k, o, 'stable');  o = [o, k(1:Kmax-numel(o))];
    end
end

function R = load_cache(file)
    R = containers.Map();
    if isfile(file), C = load(file); for i = 1:numel(C.keys), R(C.keys{i}) = C.vals{i}; end, end
end

function save_cache(file, R)
    keys = R.keys;  vals = R.values; %#ok<NASGU>
    save(file, 'keys', 'vals', '-v7.3');
end

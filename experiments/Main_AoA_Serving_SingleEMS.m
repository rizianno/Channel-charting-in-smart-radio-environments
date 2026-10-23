%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SINGLE-EMS DEPLOYMENT: AoA ESTIMATION AND EMS BEAM SELECTION WITH ONLY ONE SURFACE
%
% Same pipeline as Main_AoA_Serving_NLoS_Size, but the other EMS does not exist: its channel is
% removed from the t-SNE input AND from the served SNR (previous runs always had both EMSs).
%   - Deployment "EMS1 only" (EMS2 channel = 0) and "EMS2 only" (EMS1 channel = 0).
%   - Population (SGL_POP): fullnlos = users with no BS-UE channel; nolos = every user without LoS
%     (reflection only + full NLOS, "banana"); in both, only the users the present EMS reaches.
%     Results are reported separately for full-NLOS and reflection-only users.
%   - One 1-D t-SNE chart (azimuth at the EMS) per training acquisition, 15% anchors.
%   - Training entries (K = 1..15): nested ever-finer spread, PDF-informed (greedy coverage of the
%     anchors' angles, footprint SGL_FOOT in sin), uniform random (SGL_DRAWS draws, averaged).
%   - Per-user fusion of the K charts: avg (continuity / sum), sharp (trustworthiness - min),
%     select (the chart with the highest trustworthiness); all label-free.
%   - Beam selection: t-SNE AoA -> LUT (steering law); + R neighbours; + hierarchical search
%     warm-started at the LUT entry; searches: greedy (sequential), hierarchical, exhaustive; all
%     searches also cluster-first ([clusters]: entries within 2 beamwidths of an anchor, greedy
%     ordered by the prior density). Search sweeps up to SGL_SMAX.
%   - Metrics: AoA |error| median / RMSE / 90th percentile; SNR after steering (median, 10th
%     percentile, gain vs zero phase).
% Surfaces as in Main_AoA_Serving_NLoS_Size: 60 x 60 ray-traced, 60..240 synthesized.
% Environment: SGL_EMS ("1,2"), SGL_POP ("fullnlos,nolos"), SGL_SIZES ("60r,60,120,180,240"),
%   SGL_KS ("1:15"), SGL_RS ("2,4,6"), SGL_SMAX (30), SGL_DRAWS (10), SGL_FOOT (0.15),
%   SGL_NOISE (est_G60dB), SGL_WORKERS (48)
% Outputs: results_aoa/serving_single_ems/<noise>/ (summary.csv, aoa.csv, README.md)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();  params.target = 'aoa';
CAT   = abl_catalog();
nz    = CAT.noise(strcmp({CAT.noise.name}, envs('SGL_NOISE','est_G60dB')));
emsL  = str2double(strsplit(envs('SGL_EMS','1,2'), ','));
pops  = strtrim(strsplit(envs('SGL_POP','fullnlos,nolos'), ','));
sizes = strtrim(strsplit(envs('SGL_SIZES','60r,60,120,180,240'), ','));
Ks    = eval(envs('SGL_KS','1:15'));
Rs    = str2double(strsplit(envs('SGL_RS','2,4,6'), ','));
Smax  = str2double(envs('SGL_SMAX','30'));
nDraw = str2double(envs('SGL_DRAWS','10'));
foot  = str2double(envs('SGL_FOOT','0.15'));
nW    = str2double(envs('SGL_WORKERS','48'));
k0 = 2*pi/(3e8/30e9);  d = (3e8/30e9)/2;  lamc = 3e8/30e9;

Data = load_dataset(params);
U    = abl_users(Data, params);
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);  M0 = params.M;
Pd = squeeze(sum(sum(abs(Data.H_LOS_1bounce).^2,2),3)) / (Nb*Nr);
dd = sqrt(sum(([U.xy, params.UE_height*ones(U.n,1)] - U.BS_position).^2, 2));
los = U.has_direct & 10*log10(Pd ./ (lamc./(4*pi*dd)).^2) > -36;
rho = 10^((nz.tx_dBm - nz.noise_dBm)/10);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);  specular = -azBS;
front = U.front;
reachAll = [any(Data.H_I_RIS1 ~= 0, [2 3]), any(Data.H_I_RIS2 ~= 0, [2 3])] & front;

pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', nW); end
if ~exist(fullfile(root,'cache_ablations'),'dir'), mkdir(fullfile(root,'cache_ablations')); end
rows = {};  arows = {};
outDir = fullfile(root,'results_aoa','serving_single_ems', nz.name);
if ~exist(outDir,'dir'), mkdir(outDir); end
for pp = pops                                                   % full NLOS first: its results are written first
for e = emsL
    pop = pp{1};
    if strcmp(pop,'nolos'), cand = ~los; else, cand = ~U.has_direct; end
    usr = find(cand & reachAll(:,e));                          % population seen by this EMS
    nU  = numel(usr);  isNone = ~U.has_direct(usr);
    Hd  = reshape(permute(double(Data.H_LOS_1bounce(usr,:,:)),[2 3 1]), Nb*Nr, nU);
    xyz = [U.xy(usr,:), params.UE_height*ones(nU,1)];
    nA  = round(params.Supervision_percentage * nU);
    lab = U.az(usr, e);  perp = min(params.tsne_perplexity, round(nU/5));
    scored = (1:nU).' > nA;                                     % anchors first
    uA  = sind(lab(1:nA));
    fprintf('\n=== EMS%d only, population %s: %d users (full NLOS %d, reflection only %d), %d anchors, perplexity %d\n', ...
            e, pop, nU, nnz(isNone), nnz(~isNone), nA, perp);

    % fitted geometric model of this EMS (for the synthesized surfaces)
    Ho = double(Data.(sprintf('H_o_RIS%d',e)));  [Uo, ~, Vo] = svds(Ho, 1);
    Hi = Data.(sprintf('H_I_RIS%d',e));  Hi = Hi(usr,:,:);
    best = -inf;
    for sm = [1 -1], for sn = [1 -1], for pc = [1 -1]
        b0 = sph(elem_pos(U.EMS_position(e,:), M0, d, sm, sn), U.BS_position, k0, pc);
        cb = abs(Vo' * conj(b0)) / norm(b0);
        if cb > best, best = cb; ori = [sm sn pc]; end
    end, end, end
    P0 = elem_pos(U.EMS_position(e,:), M0, d, ori(1), ori(2));
    b0 = sph(P0, U.BS_position, k0, ori(3));
    gamma = (Uo' * Ho * conj(b0)) / (b0.' * conj(b0));
    S0 = sph(P0, xyz, k0, ori(3));  alpha = zeros(nU, Nr);
    for u = 1:nU, s0 = S0(:,u); alpha(u,:) = (s0' * squeeze(Hi(u,:,:))) / (s0' * s0); end

    for sz = sizes
        real_ = endsWith(sz{1}, 'r');  M = str2double(erase(sz{1}, 'r'));
        tag = sprintf('M%d%s', M, repmat('_raytraced', 1, real_));
        wu = 0.886 * 2 / M;  N = ceil(2/wu);  N = N + mod(N,2);
        step = ((1:N) - N/2) * 360/N;
        ub = mod(sind(specular(e)) + step/180 + 1, 2) - 1;      % beam direction of each entry (sin)
        % EMS path per entry: A(:, u, i) = (Nb*Nr) channel of user u through entry i
        if real_
            Ce = zeros(N, Nb*Nr, nU);  X = reshape(permute(Hi,[2 3 1]), M0^2, []);
            for i = 1:N
                th = exp(-1j*kron(deg2rad(step(i))*(1:M0), ones(1,M0))).';
                Ce(i,:,:) = reshape(reshape(Ho*(th.*X), Nb, Nr, nU), Nb*Nr, nU);
            end
        else
            P = elem_pos(U.EMS_position(e,:), M, d, ori(1), ori(2));
            bB = sph(P, U.BS_position, k0, ori(3));  Y = zeros(M, nU);
            for mm = 1:M, k = (mm-1)*M + (1:M); Y(mm,:) = bB(k).' * sph(P(k,:), xyz, k0, ori(3)); end
            c = exp(-1j * deg2rad(step(:)) * (1:M)) * Y;            % N x nU
            Ce = zeros(N, Nb*Nr, nU);
            for r = 1:Nr
                Ce(:, (r-1)*Nb + (1:Nb), :) = reshape(gamma*Uo, 1, Nb, 1) .* reshape(c .* alpha(:,r).', N, 1, nU);
            end
        end
        A = permute(Ce, [2 3 1]);  clear Ce
        g = zeros(nU, N);  for i = 1:N, g(:,i) = lam_max(A(:,:,i), Nb, Nr).'; end   % own EMS path gain
        snrOf = @(i) rho * lam_max(Hd + pick(A, i), Nb, Nr);
        base = snrOf(repmat(N/2, nU, 1));
        exh  = snrOf(best1(g, repmat(1:N, nU, 1)));
        ok   = scored & exh(:) > 0;
        grp  = {'full NLOS', ok & isNone; 'reflection only', ok & ~isNone};
        grp  = grp(cellfun(@nnz, grp(:,2)) > 0, :);
        hdr  = {sprintf('EMS%d only', e), pop, tag, M, N};
        add  = @(rows, train, fus, scheme, K, R, s, snr) [rows; collect1(hdr, train, fus, scheme, K, R, s, snr, base, grp)];

        % ---- prior (anchors), support, searches
        near = abs(uA - ub);
        Sup  = find(any(near <= 2*wu, 1));  Ns = numel(Sup);
        [~, o] = sort(-sum(near(:,Sup) <= wu, 1));  sq = [Sup(o), setdiff(1:N, Sup)];
        rows = add(rows, '-', '-', 'baseline (zero phase)', 0, 0, 0, base);
        rows = add(rows, '-', '-', 'exhaustive', 0, 0, N, exh);
        rows = add(rows, '-', '-', 'exhaustive[clusters]', 0, 0, Ns, snrOf(best1(g, repmat(Sup, nU, 1))));
        lutT = lut1(lab, specular(e), step);
        rows = add(rows, '-', '-', 'true AoA + LUT', 0, 0, 0, snrOf(lutT));
        for s = 1:Smax
            rows = add(rows, '-', '-', 'greedy', 0, 0, s, snrOf(best1(g, repmat(1:min(s,N), nU, 1))));
            rows = add(rows, '-', '-', 'greedy[clusters]', 0, 0, s, snrOf(best1(g, repmat(sq(1:min(s,N)), nU, 1))));
            rows = add(rows, '-', '-', 'hierarchical', 0, 0, s, snrOf(hier1(g, N, s)));
            rows = add(rows, '-', '-', 'hierarchical[clusters]', 0, 0, s, snrOf(Sup(hier1(g(:,Sup), Ns, s))));
        end

        % ---- training entries and charts (one chart per entry; cached)
        Kmax = max(Ks);
        no = nest_order(N);
        lists = {'nested', {no(1:Kmax)}; ...
                 'informed', {informed_order(uA, ub, foot, Kmax, max(2, round(0.02*nA)))}; ...
                 'random', arrayfun(@(dI) randperm(RandStream('mt19937ar','Seed',1000*dI + e), N, Kmax), 1:nDraw, 'UniformOutput', false)};
        fprintf('[%s] N = %d, support %d; nested %s; informed %s -> beams %s deg\n', tag, N, Ns, mat2str(lists{1,2}{1}), ...
                mat2str(lists{2,2}{1}), mat2str(round(asind(ub(lists{2,2}{1})))));
        ent = unique(cell2mat(cellfun(@(c) cell2mat(c), lists(:,2).', 'UniformOutput', false)));
        cacheFile = fullfile(root,'cache_ablations', sprintf('single_%s_e%d_%s_%s.mat', nz.name, e, pop, tag));
        Rmap = load_cache(cacheFile);
        need = ent(~arrayfun(@(i) isKey(Rmap, sprintf('%d', i)), ent));
        futs = parallel.FevalFuture.empty;  fk = {};  done = 0;
        for i = need
            F = logcov_features(reshape(Hd + A(:,:,i), Nb, Nr, nU), sprintf('%s|single|e%d|%s|%d', tag, e, pop, i), nz, params);
            futs(end+1) = parfeval(pool, @srv_chart_twct, 1, squareform(pdist(F)), lab, nA, perp, params.K_fraction); %#ok<SAGROW>
            fk{end+1} = sprintf('%d', i); %#ok<SAGROW>
            while numel(futs) >= 2*nW
                [j, out] = fetchNext(futs);  Rmap(fk{j}) = out;  futs(j) = [];  fk(j) = [];  done = done + 1;
            end
        end
        while ~isempty(futs)
            [j, out] = fetchNext(futs);  Rmap(fk{j}) = out;  futs(j) = [];  fk(j) = [];  done = done + 1;
        end
        if ~isempty(need), save_cache(cacheFile, Rmap); end
        fprintf('  %d charts computed, %d cached (%.1f min)\n', done, numel(ent) - done, toc(t0)/60);

        % ---- t-SNE AoA -> beam, per training list, fusion and K
        err1 = @(az) abs(mod(az - lab + 180, 360) - 180);
        for t = 1:size(lists,1)
            for K = Ks
                rr = {};  ar = {};
                for dI = 1:numel(lists{t,2})
                    li = lists{t,2}{dI}(1:K);
                    outs = Rmap(sprintf('%d', li(1)));
                    for q = 2:K, outs(q) = Rmap(sprintf('%d', li(q))); end
                    [fa, fs, fx] = fuse_variants([outs.est], [outs.tw], [outs.ct], nA, lab);
                    for fv = {'avg', fa; 'sharp', fs; 'select', fx}.'
                        az = fv{2};  lut = lut1(az, specular(e), step);
                        [~, lpos] = min(abs(ub(lut).' - ub(Sup)), [], 2);
                        rr = add(rr, lists{t,1}, fv{1}, 'tsne', K, 0, K, snrOf(lut));
                        for R = Rs
                            rr = add(rr, lists{t,1}, fv{1}, 'tsne+refine', K, R, K + R + 1, snrOf(best1(g, nbr(lut, R, N))));
                            rr = add(rr, lists{t,1}, fv{1}, 'tsne+refine[clusters]', K, R, K + R + 1, snrOf(best1(g, Sup(nbr(lpos, R, Ns)))));
                        end
                        for b = 1:(Smax - K)
                            rr = add(rr, lists{t,1}, fv{1}, 'tsne+hierarchical', K, 0, K + b, snrOf(warm1(g, lut, N, b, round(N/8))));
                            rr = add(rr, lists{t,1}, fv{1}, 'tsne+hierarchical[clusters]', K, 0, K + b, snrOf(Sup(warm1(g(:,Sup), lpos, Ns, b, max(1, round(Ns/8))))));
                        end
                        ae = err1(az);
                        for gi = 1:size(grp,1)
                            m = grp{gi,2};
                            ar(end+1,:) = [hdr, {lists{t,1}, fv{1}, K, grp{gi,1}, nnz(m), median(ae(m)), sqrt(mean(ae(m).^2)), prctile(ae(m), 90)}]; %#ok<SAGROW>
                        end
                    end
                end
                rows  = [rows;  avg_rows(rr, 14:19, [6 7 8 9 10 11 12])]; %#ok<AGROW>
                arows = [arows; avg_rows(ar, 11:13, [6 7 8 9])]; %#ok<AGROW>
            end
            q = arows(strcmp(arows(:,3), tag) & strcmp(arows(:,2), pop) & strcmp(arows(:,1), hdr{1}) & strcmp(arows(:,6), lists{t,1}) & strcmp(arows(:,7), 'select'), :);
            q = q(strcmp(q(:,9), 'full NLOS') | size(grp,1) == 1, :);
            fprintf('  %-8s select, AoA median |err| K=1..%d: %s deg\n', lists{t,1}, Kmax, mat2str(round(cell2mat(q(:,11)).', 1)));
        end
        clear A g
    end
    write_tables(outDir, rows, arows);                          % partial results after every deployment
end
end
fprintf('\n=== DONE (%.1f min) ===\n', toc(t0)/60);

%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function write_tables(outDir, rows, arows)
    T = cell2table(rows, 'VariableNames', {'deployment','population','surface','M','N','training','fusion','scheme','K','R','sweeps','group','n_users', ...
        'dsnr_median_dB','dsnr_mean_dB','dsnr_p10_dB','snr_median_dB','snr_p10_dB','se_mean'});
    TA = cell2table(arows, 'VariableNames', {'deployment','population','surface','M','N','training','fusion','K','group','n_users', ...
        'aoa_median_deg','aoa_rmse_deg','aoa_p90_deg'});
    writetable(T, fullfile(outDir,'summary.csv'));
    writetable(TA, fullfile(outDir,'aoa.csv'));
end

function r = collect1(hdr, train, fus, scheme, K, R, s, snr, base, grp)
    r = cell(size(grp,1), 19);
    for gi = 1:size(grp,1)
        m = grp{gi,2};  x = snr(m);  b = base(m);
        dd = 10*log10(x ./ b);  dd = dd(isfinite(dd));
        r(gi,:) = [hdr, {train, fus, scheme, K, R, s, grp{gi,1}, nnz(m), median(dd), mean(dd), prctile(dd,10), ...
                   median(10*log10(x)), prctile(10*log10(x(x>0)),10), mean(log2(1+x))}];
    end
end

function out = avg_rows(rr, numCols, keyCols)
% rows of several random draws -> one row per key, numeric columns averaged
    if isempty(rr), out = rr; return; end
    key = arrayfun(@(i) strjoin(cellfun(@num2str, rr(i,keyCols), 'UniformOutput', false), '|'), (1:size(rr,1)).', 'UniformOutput', false);
    [u, ia, ic] = unique(key, 'stable');
    if numel(u) == size(rr,1), out = rr; return; end
    out = rr(ia,:);  X = cell2mat(rr(:,numCols));
    for j = 1:numel(u), out(j,numCols) = num2cell(mean(X(ic==j,:), 1, 'omitnan')); end
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
% log-covariance features (same model as Main_AoA_Serving_NLoS_Size)
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

function H = pick(Ae, idx)
    [dd, n, ~] = size(Ae);
    H = reshape(Ae(:, (idx(:).'-1)*n + (1:n)), dd, n);
end

function b = best1(g, c)
% per user, the candidate entry (row of c) with the highest own-EMS gain g (n x N)
    n = size(c,1);  [~, k] = max(g((c - 1)*n + (1:n).'), [], 2);  b = c((k - 1)*n + (1:n).');
end

function c = nbr(b, R, N)
% the entry b and its R neighbours (R/2 on each side)
    c = mod(b(:) - 1 + [0, -(R/2):-1, 1:(R/2)], N) + 1;
end

function b = hier1(g, N, s)
% hierarchical search: spread over min(s,4) entries, then +-delta around the best, delta halved
    n = size(g,1);  spr = @(k) unique(round((2*(1:k) - 1) * N / (2*k)));
    b = best1(g, repmat(spr(min(s,4)), n, 1));  used = min(s,4);  delta = N/4;
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
% hierarchical search warm-started at entry b (1 sweep to measure it), +-delta, delta halved
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
% codebook entry steering toward azimuth az (steering law, lambda/2)
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
% PDF-informed entries (as in Main_AoA_Serving_NLoS_Size): greedy coverage of the anchors' angles
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

function [a, s, x] = fuse_variants(E, QT, QC, nA, lab)
% per-user fusion of K charts (as in Main_AoA_Serving_NLoS_Size)
    [n, K] = size(E);
    if K == 1, a = E; s = E; x = E; return; end
    W  = max(QC, 1e-6);  W = W ./ sum(W,2);  a = sum(W .* E, 2);
    Wt = QT - min(QT, [], 2);  z = sum(Wt,2) <= 0;  Wt(z,:) = 1;  Wt = Wt ./ sum(Wt,2);  s = sum(Wt .* E, 2);
    [~, k] = max(QT, [], 2);  x = E(sub2ind([n K], (1:n).', k));
    a(1:nA) = lab(1:nA);  s(1:nA) = lab(1:nA);  x(1:nA) = lab(1:nA);
end

function R = load_cache(file)
    R = containers.Map();
    if isfile(file), C = load(file); for i = 1:numel(C.keys), R(C.keys{i}) = C.vals{i}; end, end
end

function save_cache(file, R)
    keys = R.keys;  vals = R.values; %#ok<NASGU>
    save(file, 'keys', 'vals', '-v7.3');
end

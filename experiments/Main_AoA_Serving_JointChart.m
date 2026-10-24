%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JOINT ANGLE ESTIMATION: ONE 2-D t-SNE CHART (EMS1 AND EMS2 AZIMUTH TOGETHER), K x K TRAINING
%
% Both EMSs present. Training = real K x K grid (two nested loops): acquisition (i, j) sets EMS1 to
% codeword list1(i) and EMS2 to codeword list2(j), i, j = 1..K (K^2 acquisitions). Two lists:
%   coarse-to-fine : each EMS's first K codewords of the ever-finer order N/2, N, N/4, 3N/4, ...
%   informed       : each EMS's first K codewords of its prior-based list (greedy coverage of the
%                    anchors' angles, as Main_AoA_Serving_NLoS_Size)
% Per acquisition ONE semi-supervised t-SNE chart with 2-D output (azimuth at EMS1, at EMS2) on the
% full-NLOS users reached by at least one EMS (15% anchors; labels asind(sin(az)), continuous also
% for users behind a surface). Per user, the chart with the highest label-free trustworthiness of
% the K^2 is selected; each azimuth -> LUT -> beam of that EMS.
% Metrics (scored = not an anchor): per EMS (users it reaches): |error| of the AoA estimate and
% pointing error of the chosen beam (median, RMSE, p90), SNR served by that EMS alone; pair (users
% reached by both): sqrt(e1^2 + e2^2), SNR served by both EMSs. Schemes vs total sweeps: t-SNE
% (K^2), + R neighbours (K^2 + R + 1), + hierarchical warm start (K^2 + b), searches (greedy,
% hierarchical; full codebook and [clusters]), exhaustive. Charts cached per pair (training
% independent): a larger K later only computes the new pairs.
% Environment: SJC_SIZES ("60r,60,120,180,240"), SJC_KS ("1:8"), SJC_RS ("2,4,6"), SJC_BMAX (16,
%   hierarchical sweeps after training), SJC_FOOT (0.15), SJC_NOISE (est_G60dB), SJC_WORKERS (40)
% Outputs: results_aoa/serving_joint_chart/<noise>/summary.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();  params.target = 'aoa';
CAT   = abl_catalog();
nz    = CAT.noise(strcmp({CAT.noise.name}, envs('SJC_NOISE','est_G60dB')));
sizes = strtrim(strsplit(envs('SJC_SIZES','60r,60,120,180,240'), ','));
Ks    = eval(envs('SJC_KS','1:8'));  Kmax = max(Ks);
Rs    = str2double(strsplit(envs('SJC_RS','2,4,6'), ','));
Bmax  = str2double(envs('SJC_BMAX','16'));
foot  = str2double(envs('SJC_FOOT','0.15'));
nW    = str2double(envs('SJC_WORKERS','40'));
Smax  = Kmax^2 + max([Rs + 1, Bmax]);                           % searches up to the longest hybrid
k0 = 2*pi/(3e8/30e9);  d = (3e8/30e9)/2;

Data = load_dataset(params);
U    = abl_users(Data, params);
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);  M0 = params.M;
rho = 10^((nz.tx_dBm - nz.noise_dBm)/10);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);  specular = -azBS;

popU = find(~U.has_direct);  nU = numel(popU);                  % full NLOS (order as the earlier runs)
reach = [any(Data.H_I_RIS1(popU,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(popU,:,:) ~= 0, [2 3])] & U.front(popU,:);
Hd  = reshape(permute(double(Data.H_LOS_1bounce(popU,:,:)),[2 3 1]), Nb*Nr, nU);   % zero for full NLOS
xyz = [U.xy(popU,:), params.UE_height*ones(nU,1)];
cu  = find(any(reach,2));  nC = numel(cu);                      % chart users: reached by at least one EMS
nA  = round(params.Supervision_percentage * nC);
azTrue = U.az(popU,:);
lab = asind(sind(azTrue(cu,:)));                                % 2-D labels (front users: = azimuth)
perp = min(params.tsne_perplexity, round(nC/5));
isAnch = false(nU,1);  isAnch(cu(1:nA)) = true;
fit = fit_model(Data, U, popU, xyz, M0, d, k0);
fprintf('Joint chart: %d full-NLOS users reached by an EMS (EMS1 %d, EMS2 %d, both %d), %d anchors, perplexity %d\n', ...
        nC, nnz(reach(:,1)), nnz(reach(:,2)), nnz(all(reach,2)), nA, perp);

pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', nW); end
outDir = fullfile(root,'results_aoa','serving_joint_chart', nz.name);
if ~exist(outDir,'dir'), mkdir(outDir); end
rows = {};
for sz = sizes
    real_ = endsWith(sz{1}, 'r');  M = str2double(erase(sz{1}, 'r'));
    tag = sprintf('M%d%s', M, repmat('_raytraced', 1, real_));
    wu = 0.886 * 2 / M;  N = ceil(2/wu);  N = N + mod(N,2);
    step = ((1:N) - N/2) * 360/N;
    A = ems_paths(Data, U, popU, xyz, fit, real_, M, M0, N, step, d, k0, Nb, Nr);
    G = cellfun(@(Ae) gains(Ae, Nb, Nr), A, 'UniformOutput', false);
    lmax = @(H) lam_max(H, Nb, Nr);
    ub = arrayfun(@(e) mod(sind(specular(e)) + step/180 + 1, 2) - 1, 1:2, 'UniformOutput', false);
    no = nest_order(N);  Sup = cell(1,2);  uAe = cell(1,2);  lists = struct('c2f', {cell(1,2)}, 'informed', {cell(1,2)});
    for e = 1:2
        uA = sind(azTrue(cu(1:nA), e));  uA = uA(reach(cu(1:nA), e));  uAe{e} = uA;   % prior: anchors this EMS reaches
        Sup{e} = find(any(abs(uA - ub{e}) <= 2*wu, 1));
        lists.c2f{e} = no(1:Kmax);
        lists.informed{e} = informed_order(uA, ub{e}, foot, Kmax, max(2, round(0.02*numel(uA))));
    end
    fprintf('\n[%s] N = %d; informed EMS1 %s, EMS2 %s\n', tag, N, mat2str(lists.informed{1}), mat2str(lists.informed{2}));

    % ---- joint charts for every pair of both K x K grids (cached per pair)
    P = zeros(0,2);
    for tr = {'c2f','informed'}, [ga, gb] = ndgrid(lists.(tr{1}){1}, lists.(tr{1}){2}); P = [P; ga(:) gb(:)]; end %#ok<AGROW>
    cf = fullfile(root,'cache_ablations', sprintf('jointchart_%s_%s.mat', nz.name, tag));
    Rj = joint2d_charts(load_cache(cf), cf, unique(P,'rows','stable'), A, Hd, cu, lab, nA, perp, nz, tag, params, pool, nW, Nb, Nr, nU);
    fprintf('  charts ready (%.1f min)\n', toc(t0)/60);

    % ---- views and searches
    base = {rho * lmax(Hd + pick(A{1}, N/2*ones(nU,1))), rho * lmax(Hd + pick(A{2}, N/2*ones(nU,1)))};
    baseJ = rho * lmax(Hd + pick(A{1}, N/2*ones(nU,1)) + pick(A{2}, N/2*ones(nU,1)));
    snrE = @(e, i) rho * lmax(Hd + pick(A{e}, i));
    snrJ = @(i1, i2) rho * lmax(Hd + pick(A{1}, i1) + pick(A{2}, i2));
    exh = {best1(G{1}, repmat(1:N, nU, 1)), best1(G{2}, repmat(1:N, nU, 1))};
    mE = {reach(:,1) & ~isAnch & snrE(1, exh{1}).' > 0, reach(:,2) & ~isAnch & snrE(2, exh{2}).' > 0};
    mB = mE{1} & mE{2};                                         % pair metrics: users reached by both
    mJ = any(reach,2) & ~isAnch;
    berr = @(e, i) beam_err(ub{e}(i(:)).', azTrue(:,e), wu/2);
    hdr = {tag, M, N};
    function_rows = @(train, scheme, K, R, s, i1, i2, az) eval_rows(hdr, train, scheme, K, R, s, i1, i2, az, ...
        snrE, snrJ, base, baseJ, berr, azTrue, mE, mB, mJ);
    sq = cell(1,2);
    for e = 1:2
        [~, o] = sort(-sum(abs(uAe{e} - ub{e}(Sup{e})) <= wu, 1));  sq{e} = [Sup{e}(o), setdiff(1:N, Sup{e})];
    end
    srch = {'exhaustive',             @(e, s) exh{e}; ...
            'greedy',                 @(e, s) best1(G{e}, repmat(1:min(s,N), nU, 1)); ...
            'greedy[clusters]',       @(e, s) best1(G{e}, repmat(sq{e}(1:min(s,N)), nU, 1)); ...
            'hierarchical',           @(e, s) hier1(G{e}, N, s); ...
            'hierarchical[clusters]', @(e, s) Sup{e}(hier1(G{e}(:,Sup{e}), numel(Sup{e}), s)).'};
    for q = 1:size(srch,1)
        sw = 1:Smax;  if strcmp(srch{q,1}, 'exhaustive'), sw = 0; end
        for s = sw
            rows = [rows; function_rows('-', srch{q,1}, 0, 0, s, srch{q,2}(1, s), srch{q,2}(2, s), [])]; %#ok<AGROW>
        end
    end

    % ---- t-SNE joint chart, K x K training, per-user selection; hybrids per EMS
    for tr = {'c2f','informed'}
        for K = Ks
            [ga, gb] = ndgrid(lists.(tr{1}){1}(1:K), lists.(tr{1}){2}(1:K));  Pk = [ga(:) gb(:)];  nAcq = size(Pk,1);
            outs = Rj(sprintf('%d_%d', Pk(1,:)));
            for q = 2:nAcq, outs(q) = Rj(sprintf('%d_%d', Pk(q,:))); end
            est = select2d(outs, nA, lab);                          % nC x 2
            az = nan(nU,2);  az(cu,:) = est;
            lut = [lut1(az(:,1), specular(1), step), lut1(az(:,2), specular(2), step)];
            lpos = cell(1,2);
            for e = 1:2, [~, lpos{e}] = min(abs(ub{e}(lut(:,e)).' - ub{e}(Sup{e})), [], 2); end
            rows = [rows; function_rows(tr{1}, 'tsne', K, 0, nAcq, lut(:,1), lut(:,2), az)]; %#ok<AGROW>
            for R = Rs
                rows = [rows; function_rows(tr{1}, 'tsne+refine', K, R, nAcq + R + 1, ...
                    best1(G{1}, nbr(lut(:,1), R, N)), best1(G{2}, nbr(lut(:,2), R, N)), az)]; %#ok<AGROW>
                rows = [rows; function_rows(tr{1}, 'tsne+refine[clusters]', K, R, nAcq + R + 1, ...
                    best1(G{1}, Sup{1}(nbr(lpos{1}, R, numel(Sup{1})))), best1(G{2}, Sup{2}(nbr(lpos{2}, R, numel(Sup{2})))), az)]; %#ok<AGROW>
            end
            for b = 1:Bmax
                rows = [rows; function_rows(tr{1}, 'tsne+hierarchical', K, 0, nAcq + b, ...
                    warm1(G{1}, lut(:,1), N, b, round(N/8)), warm1(G{2}, lut(:,2), N, b, round(N/8)), az)]; %#ok<AGROW>
                w = cell(1,2);
                for e = 1:2, w{e} = Sup{e}(warm1(G{e}(:,Sup{e}), lpos{e}, numel(Sup{e}), b, max(1, round(numel(Sup{e})/8)))).'; end
                rows = [rows; function_rows(tr{1}, 'tsne+hierarchical[clusters]', K, 0, nAcq + b, w{1}, w{2}, az)]; %#ok<AGROW>
            end
            e1 = abs(az(:,1) - azTrue(:,1));  e2 = abs(az(:,2) - azTrue(:,2));
            fprintf('  %-8s K=%d (%2d acq): AoA median/RMSE EMS1 %.1f/%.1f, EMS2 %.1f/%.1f, pair %.1f/%.1f deg\n', tr{1}, K, nAcq, ...
                median(e1(mE{1})), sqrt(mean(e1(mE{1}).^2)), median(e2(mE{2})), sqrt(mean(e2(mE{2}).^2)), ...
                median(hypot(e1(mB), e2(mB))), sqrt(mean(e1(mB).^2 + e2(mB).^2)));
        end
    end
    clear A G
    write_table(outDir, rows);
end
fprintf('\n=== DONE (%.1f min) ===\n', toc(t0)/60);

%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function r = eval_rows(hdr, train, scheme, K, R, s, i1, i2, az, snrE, snrJ, base, baseJ, berr, azTrue, mE, mB, mJ)
% rows for the views EMS1, EMS2 (served by that EMS alone) and pair (served by both): SNR gain,
% beam pointing error and AoA-estimate error (median, RMSE, p90); pair errors = sqrt(e1^2 + e2^2)
    ii = {i1(:), i2(:)};  bE = {berr(1, ii{1}), berr(2, ii{2})};
    if isempty(az), aE = {nan(size(bE{1})), nan(size(bE{2}))};
    else, aE = {abs(az(:,1) - azTrue(:,1)), abs(az(:,2) - azTrue(:,2))}; end
    r = cell(3, 21);
    for e = 1:2
        r(e,:) = row1([hdr, {train, scheme, K, R, s, sprintf('EMS%d', e)}], snrE(e, ii{e}), base{e}, mE{e}, bE{e}, aE{e});
    end
    pb = hypot(bE{1}, bE{2});  pa = hypot(aE{1}, aE{2});
    rp = row1([hdr, {train, scheme, K, R, s, 'pair'}], snrJ(ii{1}, ii{2}), baseJ, mJ, nan(size(pb)), nan(size(pa)));
    rp(14:19) = {median(pb(mB)), sqrt(mean(pb(mB).^2)), prctile(pb(mB),90), median(pa(mB)), sqrt(mean(pa(mB).^2)), prctile(pa(mB),90)};
    r(3,:) = rp;
end

function r = row1(key, snr, base, m, aB, aT)
% key (9 fields), n, SNR gain median/mean, SNR median, beam error median/RMSE/p90, AoA error
% median/RMSE/p90, mean spectral efficiency, placeholder
    x = snr(m);  b = base(m);  dd = 10*log10(x(:) ./ b(:));  dd = dd(isfinite(dd));
    r = [key, {nnz(m), median(dd), mean(dd), median(10*log10(x)), ...
               median(aB(m)), sqrt(mean(aB(m).^2)), prctile(aB(m),90), median(aT(m)), sqrt(mean(aT(m).^2)), prctile(aT(m),90), ...
               mean(log2(1+x)), 0}];
end

function write_table(outDir, rows)
    T = cell2table(rows, 'VariableNames', {'surface','M','N','training','scheme','K','R','sweeps','view','n_users', ...
        'dsnr_median_dB','dsnr_mean_dB','snr_median_dB','beam_err_median_deg','beam_err_rmse_deg','beam_err_p90_deg', ...
        'tsne_err_median_deg','tsne_err_rmse_deg','tsne_err_p90_deg','se_mean','x'});
    T.x = [];
    writetable(T, fullfile(outDir,'summary.csv'));
end

function est = select2d(outs, nA, lab)
% per user, the 2-D estimate of the chart with the highest label-free trustworthiness
    E = cat(3, outs.est);  QT = [outs.tw];  [n, ~] = size(QT);
    [~, k] = max(QT, [], 2);
    est = zeros(n, 2);
    for c = 1:2, Ec = squeeze(E(:,c,:));  if n == 1, Ec = Ec.'; end, est(:,c) = Ec(sub2ind(size(Ec), (1:n).', k)); end
    est(1:nA,:) = lab(1:nA,:);
end

function Rmap = joint2d_charts(Rmap, cacheFile, P, A, Hd, cu, lab, nA, perp, nz, tag, params, pool, nW, Nb, Nr, nU)
% one 2-D chart per pair (EMS1 codeword, EMS2 codeword): both azimuths as anchor labels
    futs = parallel.FevalFuture.empty;  fk = {};  nCh = 0;
    for q = 1:size(P,1)
        key = sprintf('%d_%d', P(q,:));
        if isKey(Rmap, key), continue; end
        H = reshape(Hd + A{1}(:,:,P(q,1)) + A{2}(:,:,P(q,2)), Nb, Nr, nU);
        if P(q,1) == P(q,2), skey = sprintf('%s|%d', tag, P(q,1)); else, skey = sprintf('%s|%d_%d', tag, P(q,:)); end
        F = logcov_features(H, skey, nz, params);
        futs(end+1) = parfeval(pool, @srv_chart_twct, 1, squareform(pdist(F(cu,:))), lab, nA, perp, params.K_fraction); %#ok<AGROW>
        fk{end+1} = key; %#ok<AGROW>
        while numel(futs) >= 2*nW, [j, out] = fetchNext(futs); Rmap(fk{j}) = out; futs(j) = []; fk(j) = []; nCh = nCh + 1; end
    end
    while ~isempty(futs)
        [j, out] = fetchNext(futs);  Rmap(fk{j}) = out;  futs(j) = [];  fk(j) = [];  nCh = nCh + 1;
        if mod(nCh, 50) == 0, save_cache(cacheFile, Rmap); end
    end
    if nCh > 0, save_cache(cacheFile, Rmap); end
    fprintf('  joint 2-D charts: %d computed, %d cached\n', nCh, size(P,1) - nCh);
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

function err = beam_err(u, az, hw)
% |pointing error| [deg] of beams at sin = u (in [-1, 1)); at lambda/2 the copy u -+ 2 is radiated
% only if within hw (half a beamwidth) of the visible region (edge beams), as Main_AoA_Serving_Deployments
    ut = sind(az);  C = [u, u - 2, u + 2];  C(abs(C) > 1 + hw) = NaN;
    [~, k] = min(abs(C - ut), [], 2);  u2 = C(sub2ind(size(C), (1:numel(u)).', k));
    err = abs(asind(min(max(u2, -1), 1)) - az);
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
% incremental hierarchical search (as Main_AoA_Serving_Deployments): first probes = first entries
% of the coarse-to-fine order, then +-delta around the best so far, delta halved
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
% coarse-to-fine order: N/2, N, N/4, 3N/4, N/8, ...
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

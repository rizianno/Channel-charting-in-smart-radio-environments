%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FULL-NLOS USERS: EMS BEAM SELECTION vs NUMBER OF SWEEPS, FOR SURFACE SIZES 60..240
%
% Same comparison as Main_AoA_Serving_Sweeps (sequential / spread / hierarchical sweeping,
% spread + refinement, t-SNE AoA + LUT, t-SNE + refinement, exhaustive; all vs number of
% sweeps, codebook fixed), restricted to the users in FULL NLOS: no BS-UE channel at all
% (no LoS and no reflection), so the EMSs are their only link.
%   - t-SNE charts are built from these users only: one 1-D chart per EMS with the full-NLOS
%     users that EMS reaches (in front of it, non-zero incident channel), 15% of them as
%     anchors, perplexity min(150, n/5) (small sets; larger perplexity makes t-SNE diverge).
%   - Surface sizes M x M (lambda/2 spacing, steering along the m axis only, the same phase
%     along the vertical n axis): 60 = the ray-traced data; 60, 120, 180, 240 synthesized:
%       incident channel  h(k) = alpha_u,r * exp(-j k0 |p_k - p_UE|) / |p_k - p_UE|
%       EMS -> BS channel H_o = gamma * a_BS * b^T, b(k) = exp(-j k0 |p_k - p_BS|) / |p_k - p_BS|
%     with alpha (per user and UE antenna), a_BS and gamma fitted on the 60 x 60 ray-traced
%     channels (validated in the README: correlation and best-beam agreement). Larger surfaces
%     get their gain and narrower beams from the geometry (near field included); partial
%     blockage of the larger aperture is not modelled.
%   - Codebook N = ceil(2 / (0.886 * 2/M)) (even): -3 dB beam crossings over sin(theta) in
%     [-1, 1]; M = 60/120/180/240 -> N = 68/136/204/272, phase steps over -180..180.
% Selection by measurement uses each EMS path's own gain (spatially separated at the BS),
% noise-free; t-SNE uses estimated channels (SRV_NOISE); the served SNR uses the true channel.
%
% SRV_POP = nolos ("banana experiment"): the same, but the population (charts, anchors, scores)
% is every user WITHOUT LoS: reflection only (BS link only via a wall; direct-channel power
% more than 36 dB below free space) plus full NLOS. Their BS channel is then part of both the
% t-SNE input and the served SNR. Outputs in results_aoa/serving_banana/<noise>/.
% SRV_TRAIN = same (default): training sweep k sets both EMSs to training entry k (K sweeps);
%           = grid: both EMSs stepped independently over the K training entries (K^2 sweeps).
% Environment: SRV_TRAIN (same | grid), SRV_POP (fullnlos | nolos), SRV_NOISE (est_G60dB), SRV_SIZES ("60r,60,120,180,240"; "60r" = ray-traced),
%   SRV_KS ("1:8"), SRV_RS ("2,4,6"), SRV_SMAX (40), SRV_WORKERS (48)
% Outputs: results_aoa/serving_nlos_size/<noise>/ (README.md, summary.csv, validation.csv, figures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();  params.target = 'aoa';
CAT  = abl_catalog();
nz   = CAT.noise(strcmp({CAT.noise.name}, envs('SRV_NOISE','est_G60dB')));
sizes = strtrim(strsplit(envs('SRV_SIZES','60r,60,120,180,240'), ','));
Ks   = eval(envs('SRV_KS','1:8'));
Rs   = str2double(strsplit(envs('SRV_RS','2,4,6'), ','));
Smax = str2double(envs('SRV_SMAX','40'));
nW   = str2double(envs('SRV_WORKERS','48'));
k0   = 2*pi/(3e8/30e9);  d = (3e8/30e9)/2;

%% ------------------------------------------------------------------------
% DATA, FULL-NLOS USERS, FIT OF THE GEOMETRIC MODEL ON THE 60 x 60 DATA
%% ------------------------------------------------------------------------
Data = load_dataset(params);
U    = abl_users(Data, params);
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);  M0 = params.M;
pop  = envs('SRV_POP','fullnlos');
train = envs('SRV_TRAIN','same');
lamc = 3e8/30e9;
Pd   = squeeze(sum(sum(abs(Data.H_LOS_1bounce).^2,2),3)) / (Nb*Nr);
dd   = sqrt(sum(([U.xy, params.UE_height*ones(U.n,1)] - U.BS_position).^2, 2));
los  = U.has_direct & 10*log10(Pd ./ (lamc./(4*pi*dd)).^2) > -36;
if strcmp(pop,'nolos'), nlos = find(~los); else, nlos = find(~U.has_direct); end   % population
nU   = numel(nlos);
isNone = ~U.has_direct(nlos);                                % full NLOS within the population
Hd   = reshape(permute(double(Data.H_LOS_1bounce(nlos,:,:)),[2 3 1]), Nb*Nr, nU);   % zero for full NLOS
xyz  = [U.xy(nlos,:), params.UE_height*ones(nU,1)];
rho  = 10^((nz.tx_dBm - nz.noise_dBm)/10);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);
specular = -azBS;
fprintf('Population %s: %d users (full NLOS %d, reflection only %d)\n', pop, nU, nnz(isNone), nnz(~isNone));

% orientation of the element grid (m axis along +-y, n axis along +-z) and fitted model
fit = struct();  val = {};
for e = 1:2
    Ho = double(Data.(sprintf('H_o_RIS%d',e)));
    [Uo, So, Vo] = svds(Ho, 1);
    Hi = Data.(sprintf('H_I_RIS%d',e));  Hi = Hi(nlos,:,:);       % nU x M0^2 x Nr
    best = -inf;                                                    % grid orientation and phase sign
    for sm = [1 -1], for sn = [1 -1], for pc = [1 -1]
        P0 = elem_pos(U.EMS_position(e,:), M0, d, sm, sn);
        b0 = sph(P0, U.BS_position, k0, pc);
        cb = abs(Vo' * conj(b0)) / norm(b0);                        % Ho = u s v' ~ gamma a b.'  ->  v ~ conj(b)
        if cb > best, best = cb; ori = [sm sn pc]; end
    end, end, end
    P0 = elem_pos(U.EMS_position(e,:), M0, d, ori(1), ori(2));
    b0 = sph(P0, U.BS_position, k0, ori(3));
    gamma = (Uo' * Ho * conj(b0)) / (b0.' * conj(b0));             % Ho ~ gamma * a * b.'  (a = Uo)
    S0 = sph(P0, xyz, k0, ori(3));                                 % M0^2 x nU
    alpha = zeros(nU, Nr);  cor = zeros(nU,1);
    for u = 1:nU
        h = squeeze(Hi(u,:,:));                                    % M0^2 x Nr
        s = S0(:,u);
        alpha(u,:) = (s' * h) / (s' * s);
        if any(h(:)), cor(u) = norm(s' * h) / (norm(s) * norm(h)); end
    end
    fit(e).ori = ori;  fit(e).a = Uo;  fit(e).gamma = gamma;  fit(e).alpha = alpha;
    fit(e).cor = cor;  fit(e).corBS = best;
    reach = any(alpha ~= 0, 2) & U.front(nlos,e);
    val(end+1,:) = {e, best, median(cor(reach)), 100*mean(cor(reach) > 0.9), nnz(reach)}; %#ok<SAGROW>
    fprintf('EMS%d: orientation m %+d, n %+d, phase sign %+d; |corr| EMS->BS %.3f; UE->EMS median %.3f, >0.9 for %.0f%% of %d users\n', ...
            e, ori(1), ori(2), ori(3), best, median(cor(reach)), 100*mean(cor(reach) > 0.9), nnz(reach));
end

%% ------------------------------------------------------------------------
% CHART STRUCTURE (per EMS, full-NLOS users it reaches)
%% ------------------------------------------------------------------------
reachE = [any(Data.H_I_RIS1(nlos,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(nlos,:,:) ~= 0, [2 3])] & U.front(nlos,:);
charts = struct('users',{}, 'nA',{}, 'labels',{}, 'perp',{});
isAnch = false(nU,2);
for e = 1:2
    m = find(reachE(:,e));
    nA = round(params.Supervision_percentage * numel(m));
    charts(e).users  = m;                                    % indices into nlos (anchors first)
    charts(e).nA     = nA;
    charts(e).labels = U.az(nlos(m), e);
    charts(e).perp   = min(params.tsne_perplexity, round(numel(m)/5));
    isAnch(m(1:nA), e) = true;
end
scored = any(reachE,2) & ~any(isAnch,2);                    % not an anchor in either chart
fprintf('Chart EMS1: %d users (%d anchors, perplexity %d); EMS2: %d (%d, %d); scored: %d\n', ...
    numel(charts(1).users), charts(1).nA, charts(1).perp, numel(charts(2).users), charts(2).nA, charts(2).perp, nnz(scored));

pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', nW); end
if ~exist(fullfile(root,'cache_ablations'),'dir'), mkdir(fullfile(root,'cache_ablations')); end
rows = {};  vrows = {};
for sz = sizes
    real_ = endsWith(sz{1}, 'r');  M = str2double(erase(sz{1}, 'r'));
    tag = sprintf('M%d%s', M, repmat('_raytraced', 1, real_));
    wu  = 0.886 * 2 / M;  N = ceil(2/wu);  N = N + mod(N,2);
    step = ((1:N) - N/2) * 360/N;
    spread = @(s) unique(round((2*(1:s) - 1) * N / (2*s)));
    fprintf('\n[%s] beamwidth (sin) %.4f -> N = %d\n', tag, wu, N);

    % per-entry EMS paths: c{e}(i,u) = scalar EMS response; path = gamma * a * c * alpha.'
    C = cell(1,2);
    for e = 1:2
        if real_
            Hi = Data.(sprintf('H_I_RIS%d',e));  Hi = Hi(nlos,:,:);
            Ho = double(Data.(sprintf('H_o_RIS%d',e)));
            X  = reshape(permute(Hi,[2 3 1]), M0^2, []);
            C{e} = zeros(N, Nb*Nr, nU);
            for i = 1:N
                th = exp(-1j*kron(deg2rad(step(i))*(1:M0), ones(1,M0))).';
                C{e}(i,:,:) = reshape(reshape(Ho*(th.*X), Nb, Nr, nU), Nb*Nr, nU);
            end
        else
            P  = elem_pos(U.EMS_position(e,:), M, d, fit(e).ori(1), fit(e).ori(2));
            b  = sph(P, U.BS_position, k0, fit(e).ori(3));
            Y  = zeros(M, nU);                                      % sum over the column n
            for mm = 1:M
                k = (mm-1)*M + (1:M);
                Y(mm,:) = b(k).' * sph(P(k,:), xyz, k0, fit(e).ori(3));
            end
            TH = exp(-1j * deg2rad(step(:)) * (1:M));               % N x M
            c  = TH * Y;                                            % N x nU
            path = fit(e).gamma * reshape(fit(e).a, Nb, 1, 1);      % Nb x 1
            C{e} = zeros(N, Nb*Nr, nU);
            for r = 1:Nr
                C{e}(:, (r-1)*Nb + (1:Nb), :) = reshape(path, 1, Nb, 1) .* reshape(c .* fit(e).alpha(:,r).', N, 1, nU);
            end
        end
    end
    A = cellfun(@(Ce) permute(Ce, [2 3 1]), C, 'UniformOutput', false);   % (Nb*Nr) x nU x N
    clear C
    Gown = cellfun(@(Ae) gain_all(Ae, Nb, Nr), A, 'UniformOutput', false);
    lmax = @(H) lam_max(H, Nb, Nr);
    snrOf = @(i1, i2) rho * lmax(Hd + pick(A{1}, i1) + pick(A{2}, i2));
    base  = snrOf(repmat(N/2,nU,1), repmat(N/2,nU,1));
    [e1, e2] = deal(repmat(1:N,nU,1));  [i1, i2] = best_of(Gown, {e1, e2});
    exh = snrOf(i1, i2);
    ok = scored & (exh(:) > 0);
    grp = {'full NLOS, reachable', ok & isNone; 'reachable by both EMSs', ok & isNone & all(reachE,2); ...
           'reachable by one EMS', ok & isNone & xor(reachE(:,1), reachE(:,2))};
    if strcmp(pop,'nolos')
        grp = [grp; {'reflection only, reachable', ok & ~isNone; 'no LoS (all), reachable', ok}];
    end
    sel = @(cand) best_of(Gown, cand);
    add = @(rows, scheme, R, s, snr) [rows; collect(tag, M, N, scheme, R, s, snr, base, grp)];
    rows = add(rows, 'baseline', 0, 0, base);
    rows = add(rows, 'exhaustive', 0, N, exh);
    for s = 1:Smax                                                   % beyond N sweeps: = exhaustive
        [i1, i2] = sel({repmat(1:min(s,N),nU,1), repmat(1:min(s,N),nU,1)});        rows = add(rows, 'sequential', 0, s, snrOf(i1,i2));
        sp = spread(min(s,N));  [i1, i2] = sel({repmat(sp,nU,1), repmat(sp,nU,1)});  rows = add(rows, 'spread', 0, s, snrOf(i1,i2));
        [i1, i2] = hierarchical(Gown, N, s, spread);                 rows = add(rows, 'hierarchical', 0, s, snrOf(i1,i2));
        for R = Rs
            if s - R < 1, continue; end
            sp = spread(min(s - R, N));  [b1, b2] = sel({repmat(sp,nU,1), repmat(sp,nU,1)});
            [i1, i2] = sel({nbr(b1, R, N), nbr(b2, R, N)});          rows = add(rows, 'spread+refine', R, s, snrOf(i1,i2));
        end
    end
    lutT = lut_physics(U.az(nlos,:), specular, step);
    rows = add(rows, 'true AoA + LUT', 0, 0, snrOf(lutT(:,1), lutT(:,2)));

    % ---- t-SNE training charts: per EMS, features of the full-NLOS users
    ctag = tag;  if strcmp(pop,'nolos'), ctag = ['nolos_' tag]; end
    cacheFile = fullfile(root,'cache_ablations', sprintf('nlos_size_%s_%s.mat', nz.name, ctag));
    Rmap = load_cache(cacheFile);
    trainPairs = @(K) pairs_of(spread(K), train);              % [entry EMS1, entry EMS2] per sweep
    Ptr  = unique(cell2mat(arrayfun(@(K) trainPairs(K), Ks, 'UniformOutput', false).'), 'rows', 'stable');
    todo = {};
    for q = 1:size(Ptr,1)
        for e = 1:2
            key = pair_key(nz.name, ctag, e, Ptr(q,:));
            if ~isKey(Rmap, key), todo(end+1,:) = {key, Ptr(q,:), e}; end %#ok<SAGROW>
        end
    end
    fprintf('  %d training configurations (%s), %d charts to compute\n', size(Ptr,1), train, size(todo,1));
    Fcache = containers.Map();
    futs = parallel.FevalFuture.empty;  fk = {};  nx = 1;  done = 0;
    while nx <= size(todo,1) || ~isempty(futs)
        while nx <= size(todo,1) && numel(futs) < 2*nW
            [key, pr, e] = todo{nx,:};  nx = nx + 1;
            fkey = sprintf('%d_%d', pr);
            if ~isKey(Fcache, fkey)
                H = Hd + A{1}(:,:,pr(1)) + A{2}(:,:,pr(2));          % EMS1 on entry pr(1), EMS2 on pr(2)
                if pr(1) == pr(2), skey = sprintf('%s|%d', ctag, pr(1)); else, skey = sprintf('%s|%d_%d', ctag, pr); end
                Fcache(fkey) = logcov_features(reshape(H, Nb, Nr, nU), skey, nz, params);
            end
            F = Fcache(fkey);  ch = charts(e);
            D = squareform(pdist(F(ch.users,:)));
            futs(end+1) = parfeval(pool, @abl_chart, 1, D, ch.labels, ch.nA, true(numel(ch.users),1), ...
                                   ch.perp, params.K_fraction); %#ok<SAGROW>
            fk{end+1} = key; %#ok<SAGROW>
        end
        [j, out] = fetchNext(futs);
        Rmap(fk{j}) = out;  futs(j) = [];  fk(j) = [];  done = done + 1;
        if mod(done,20) == 0 || (nx > size(todo,1) && isempty(futs))
            save_cache(cacheFile, Rmap);
            fprintf('  %d/%d charts (%.1f min)\n', done, size(todo,1), toc(t0)/60);
        end
    end
    clear Fcache

    % ---- t-SNE AoA -> LUT (-> refinement)
    for K = Ks
        Pk = trainPairs(K);  azHat = nan(nU,2);  diverged = 0;
        for e = 1:2
            outs = Rmap(pair_key(nz.name, ctag, e, Pk(1,:)));
            for s = 2:size(Pk,1), outs(s) = Rmap(pair_key(nz.name, ctag, e, Pk(s,:))); end
            E3 = cat(3, outs.est);  ch = charts(e);
            if size(Pk,1) == 1, est = E3; else, [~, est] = abl_fuse(E3, [outs.q_free], ch.nA, ch.labels); end
            diverged = max(diverged, mean(abs(est) > 180));
            azHat(ch.users, e) = est;
        end
        lut = lut_physics(azHat, specular, step);
        nK  = size(Pk,1);                                            % training sweeps used
        rows = add(rows, 'tsne', 0, nK, snrOf(lut(:,1), lut(:,2)));
        for R = Rs
            [i1, i2] = sel({nbr(lut(:,1), R, N, true), nbr(lut(:,2), R, N, true)});
            rows = add(rows, 'tsne+refine', R, nK + R + 1, snrOf(i1, i2));
        end
        % t-SNE warm-started hierarchical search: measure the LUT entry, then test +-delta around
        % the current best and halve delta (delta0 = N/8); column R holds K for this scheme
        for b = 1:(Smax - nK)
            [i1, i2] = warm_hier(Gown, lut, N, b, round(N/8));
            rows = add(rows, 'tsne+hierarchical', nK, nK + b, snrOf(i1, i2));
        end
        err = mod(azHat - U.az(nlos,:) + 180, 360) - 180;
        r1 = reachE(:,1) & ~isAnch(:,1);  r2 = reachE(:,2) & ~isAnch(:,2);
        offb = [mean(cdist(lut(r1,1), lutT(r1,1), N) > 1), mean(cdist(lut(r2,2), lutT(r2,2), N) > 1)];
        vrows(end+1,:) = {tag, sprintf('K=%d (%d sweeps): AoA median |err| %.1f/%.1f deg, RMSE %.1f/%.1f; LUT >1 beam off %.0f%%/%.0f%%; diverged %.1f%%', ...
            K, nK, median(abs(err(r1,1))), median(abs(err(r2,2))), sqrt(mean(err(r1,1).^2)), sqrt(mean(err(r2,2).^2)), ...
            100*offb(1), 100*offb(2), 100*diverged)}; %#ok<SAGROW>
        fprintf('  %s\n', vrows{end,2});
    end
    if ~real_ && M == M0, vrows = [vrows; beam_agreement(tag, Data, nlos, reachE, Gown, step, M0, Nb, Nr)]; end %#ok<AGROW>
    clear A Gown
end

T = cell2table(rows, 'VariableNames', {'surface','M','N','scheme','R','sweeps','group','n_users', ...
    'dsnr_median_dB','dsnr_mean_dB','dsnr_p10_dB','snr_median_dB','snr_p10_dB','se_mean','dse_mean'});
fold = 'serving_nlos_size';  if strcmp(pop,'nolos'), fold = 'serving_banana'; end
if strcmp(train,'grid'), fold = [fold '_gridTraining']; end
outDir = fullfile(root,'results_aoa', fold, nz.name);
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir,'summary.csv'));
V = cell2table(vrows, 'VariableNames', {'surface','note'});  writetable(V, fullfile(outDir,'validation_and_aoa.csv'));
for sz = unique(T.surface,'stable').'
    plot_sweeps(T(strcmp(T.surface,sz{1}),:), sz{1}, Rs, outDir, 'full NLOS, reachable');
    if strcmp(pop,'nolos')
        plot_sweeps(T(strcmp(T.surface,sz{1}),:), sz{1}, Rs, outDir, 'reflection only, reachable');
    end
end
plot_sizes(T, outDir);
write_readme(outDir, T, V, val, fit, charts, nU, nz, Ks, Rs, toc(t0)/60, pop, nnz(isNone), train);
fprintf('\n=== DONE (%.1f min) ===\n', toc(t0)/60);


%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function P = elem_pos(c, M, d, sm, sn)
% element positions (M^2 x 3), index k = (m-1) M + n: m along +-y (steering), n along +-z
    o = ((1:M) - (M+1)/2) * d;
    [n_, m_] = meshgrid(o, o);                  % rows: m, cols: n
    m_ = m_.';  n_ = n_.';                      % column-major over n first -> k = (m-1)M + n
    P = c + [zeros(M^2,1), sm*m_(:), sn*n_(:)];
end

function S = sph(P, Q, k0, pc)
% spherical-wave response exp(-j pc k0 r)/r between element positions P (K x 3) and points Q
% (q x 3); pc = +-1 is the phase-sign convention of the data (fitted)
    r = sqrt((P(:,1) - Q(:,1).').^2 + (P(:,2) - Q(:,2).').^2 + (P(:,3) - Q(:,3).').^2);
    S = exp(-1j*pc*k0*r) ./ r;
end

function F = logcov_features(H, key, noise, params)
% log-covariance features of every user from its channel H (Nb x Nr x U), same model as
% abl_state_features (estimated channels with pilot gain, or noiseless with floor)
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

function G = gain_all(Ae, Nb, Nr)
    G = zeros(size(Ae,2), size(Ae,3));
    for i = 1:size(Ae,3), G(:,i) = lam_max(Ae(:,:,i), Nb, Nr).'; end
end

function H = pick(Ae, idx)
    [dd, n, ~] = size(Ae);
    H = reshape(Ae(:, (idx(:).'-1)*n + (1:n)), dd, n);
end

function [i1, i2] = best_of(Gown, cand)
    out = cell(1,2);
    for e = 1:2
        c = cand{e};  n = size(c,1);
        g = Gown{e}((c - 1)*n + (1:n).');
        [~, k] = max(g, [], 2);
        out{e} = c((k - 1)*n + (1:n).');
    end
    [i1, i2] = out{:};
end

function c = nbr(b, R, N, withCentre)
    if nargin < 4, withCentre = false; end
    off = [-(R/2):-1, 1:(R/2)];  if withCentre, off = [0 off]; end
    c = mod(b(:) - 1 + off, N) + 1;
end

function dd = cdist(a, b, N)
    dd = abs(a - b);  dd = min(dd, N - dd);
end

function [i1, i2] = hierarchical(Gown, N, s, spread)
    n = size(Gown{1},1);  out = cell(1,2);
    for e = 1:2
        g = Gown{e};
        if s <= 4
            c = repmat(spread(s), n, 1);  [~,k] = max(g((c-1)*n + (1:n).'), [], 2);  out{e} = c((k-1)*n + (1:n).');
            continue;
        end
        c = repmat(spread(4), n, 1);  [~,k] = max(g((c-1)*n + (1:n).'), [], 2);  b = c((k-1)*n + (1:n).');
        used = 4;  delta = N/4;
        while used < s && delta > 1
            delta = ceil(delta/2);
            for sgn = [-1 1]
                if used >= s, break; end
                cand = mod(b - 1 + sgn*delta, N) + 1;
                better = g((cand-1)*n + (1:n).') > g((b-1)*n + (1:n).');
                b(better) = cand(better);  used = used + 1;
            end
        end
        out{e} = b;
    end
    [i1, i2] = out{:};
end

function P = pairs_of(tr, mode)
% training sweeps: [entry of EMS1, entry of EMS2] per sweep
    if strcmp(mode,'grid'), [a, b] = ndgrid(tr, tr); P = [a(:) b(:)];
    else, P = [tr(:) tr(:)]; end
end

function k = pair_key(nzName, ctag, e, pr)
% cache key of one chart; same-entry sweeps keep the key of the earlier runs (cache reused)
    if pr(1) == pr(2), k = sprintf('%s|%s|e%d|%d', nzName, ctag, e, pr(1));
    else, k = sprintf('%s|%s|e%d|%d_%d', nzName, ctag, e, pr(1), pr(2)); end
end

function t = trtxt(train)
    if strcmp(train,'grid'), t = 'K x K grid: both EMSs stepped independently, K^2 training sweeps';
    else, t = 'sweep k sets both EMSs to training entry k, K training sweeps'; end
end

function [i1, i2] = warm_hier(Gown, lut, N, budget, delta0)
% hierarchical search per EMS starting from the LUT entry: 1 sweep to measure it, then the two
% entries at +-delta around the current best (one sweep each), delta halved down to 1
    n = size(Gown{1},1);  out = cell(1,2);
    for e = 1:2
        g = Gown{e};  b = lut(:,e);  used = 1;  delta = delta0;
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
        out{e} = b;
    end
    [i1, i2] = out{:};
end

function sel = lut_physics(az, specular, step)
    sel = zeros(size(az,1), 2);
    for e = 1:2
        want = mod(180*(sind(az(:,e)) - sind(specular(e))) + 180, 360) - 180;
        [~, sel(:,e)] = min(abs(mod(want - step + 180, 360) - 180), [], 2);
        sel(isnan(az(:,e)), e) = numel(step)/2;
    end
end

function r = collect(tag, M, N, scheme, R, s, snr, base, grp)
    r = cell(size(grp,1), 15);
    for g = 1:size(grp,1)
        m = grp{g,2};  x = snr(m);  b = base(m);
        dd = 10*log10(x ./ b);  dd = dd(isfinite(dd));
        r(g,:) = {tag, M, N, scheme, R, s, grp{g,1}, nnz(m), median(dd), mean(dd), prctile(dd,10), ...
                  median(10*log10(x)), prctile(10*log10(x(x>0)),10), mean(log2(1+x)), mean(log2(1+x) - log2(1+b))};
    end
end

function rows = beam_agreement(tag, Data, nlos, reachE, GownSyn, step, M0, Nb, Nr)
% validation of the synthesized 60 x 60 channels against the ray-traced ones: does the best
% entry agree, and how far is its gain from the ray-traced best gain
    rows = {};
    for e = 1:2
        Hi = Data.(sprintf('H_I_RIS%d',e));  Hi = Hi(nlos,:,:);
        Ho = double(Data.(sprintf('H_o_RIS%d',e)));
        X  = reshape(permute(Hi,[2 3 1]), M0^2, []);
        G  = zeros(numel(nlos), numel(step));
        for i = 1:numel(step)
            th = exp(-1j*kron(deg2rad(step(i))*(1:M0), ones(1,M0))).';
            G(:,i) = lam_max(Ho*(th.*X), Nb, Nr).';
        end
        r = reachE(:,e);
        [gr, ir] = max(G(r,:), [], 2);  [gs, is] = max(GownSyn{e}(r,:), [], 2);
        gsr = G(r,:);  gAtSyn = gsr(sub2ind(size(gsr), (1:nnz(r)).', is));
        rows(end+1,:) = {tag, sprintf(['validation EMS%d (%d users): best entry same as ray-traced for %.0f%%, within 1 entry for %.0f%%; ' ...
            'ray-traced gain at the synthetic best entry: median %.1f dB below the ray-traced best; best gain synthetic vs ray-traced: median %+.1f dB'], ...
            e, nnz(r), 100*mean(is == ir), 100*mean(min(abs(is-ir), numel(step)-abs(is-ir)) <= 1), ...
            median(10*log10(gr ./ gAtSyn)), median(10*log10(gs ./ gr)))}; %#ok<AGROW>
        fprintf('  %s\n', rows{end,2});
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

function plot_sweeps(T, tag, Rs, outDir, grp)
    met = {'dsnr_median_dB','median SNR gain vs zero phase [dB]'; 'dsnr_p10_dB','10th percentile SNR gain [dB] (worst 10%)'; ...
           'snr_median_dB','median SNR [dB]'};
    blues = [0.30 0.75 0.93; 0.00 0.45 0.74; 0.00 0.20 0.45];  orng = [1.00 0.75 0.40; 0.95 0.55 0.15; 0.70 0.30 0.00];
    f = figure('Position',[50 50 1800 480]);
    tl = tiledlayout(f, 1, size(met,1), 'TileSpacing','compact');
    for k = 1:size(met,1)
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        q = @(sc, R) sortrows(T(strcmp(T.scheme,sc) & T.R==R & strcmp(T.group,grp),:), 'sweeps');
        r = q('sequential',0);   plot(ax, r.sweeps, r.(met{k,1}), ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5, 'DisplayName','sequential sweep');
        r = q('spread',0);       plot(ax, r.sweeps, r.(met{k,1}), '-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.8, 'DisplayName','spread sweep');
        r = q('hierarchical',0); plot(ax, r.sweeps, r.(met{k,1}), '-', 'Color',[0.49 0.18 0.56], 'LineWidth',1.8, 'DisplayName','hierarchical search');
        for i = 1:numel(Rs)
            r = q('spread+refine',Rs(i)); plot(ax, r.sweeps, r.(met{k,1}), '--', 'Color',orng(i,:), 'LineWidth',1.2, 'DisplayName', sprintf('spread sweep + %d neighbours', Rs(i)));
        end
        r = q('tsne',0); plot(ax, r.sweeps, r.(met{k,1}), '-o', 'Color',[0 0.45 0.74], 'LineWidth',2, 'MarkerSize',5, 'MarkerFaceColor',[0 0.45 0.74], 'DisplayName','t-SNE AoA + LUT');
        for i = 1:numel(Rs)
            r = q('tsne+refine',Rs(i)); plot(ax, r.sweeps, r.(met{k,1}), '-s', 'Color',blues(i,:), 'LineWidth',1.5, 'MarkerSize',5, 'DisplayName', sprintf('t-SNE AoA + LUT + %d neighbours', Rs(i)));
        end
        r = q('true AoA + LUT',0);  yline(ax, r.(met{k,1}), ':', 'Color',[0 0.45 0.74], 'LineWidth',1.5, 'DisplayName','true AoA + LUT');
        r = q('exhaustive',0);      yline(ax, r.(met{k,1}), '--k', 'LineWidth',1.5, 'DisplayName', sprintf('exhaustive (%d sweeps)', r.sweeps));
        r = q('baseline',0);        yline(ax, r.(met{k,1}), '-', 'Color',[0.6 0.6 0.6], 'DisplayName','zero phase');
        xlim(ax, [0 40]);  xlabel(ax,'number of sweeps');  ylabel(ax, met{k,2});
        if k == size(met,1), legend(ax,'Location','eastoutside','FontSize',8); end
    end
    title(tl, sprintf('%s users, surface %s (codebook %d entries)', grp, strrep(tag,'_','\_'), T.N(1)));
    fn = sprintf('sweeps_%s.png', tag);
    if ~startsWith(grp,'full'), fn = strrep(sprintf('sweeps_%s_%s.png', tag, strtok(grp,',')), ' ', '_'); end
    exportgraphics(f, fullfile(outDir, fn), 'Resolution',110);  close(f);
end

function plot_sizes(T, outDir)
% per surface size: median and 10th-percentile SNR gain vs sweeps, hierarchical vs t-SNE (+6)
    grp = 'full NLOS, reachable';  tags = unique(T.surface,'stable');  cols = lines(numel(tags));
    met = {'dsnr_median_dB','median SNR gain [dB]'; 'dsnr_p10_dB','10th percentile SNR gain [dB]'; 'snr_median_dB','median SNR [dB]'};
    f = figure('Position',[50 50 1800 480]);
    tl = tiledlayout(f, 1, size(met,1), 'TileSpacing','compact');
    for k = 1:size(met,1)
        ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
        for t = 1:numel(tags)
            q = @(sc, R) sortrows(T(strcmp(T.surface,tags{t}) & strcmp(T.scheme,sc) & T.R==R & strcmp(T.group,grp),:), 'sweeps');
            r = q('hierarchical',0); plot(ax, r.sweeps, r.(met{k,1}), '-', 'Color',cols(t,:), 'LineWidth',1.6, 'DisplayName', sprintf('%s hierarchical', strrep(tags{t},'_',' ')));
            r = q('tsne+refine',6);  plot(ax, r.sweeps, r.(met{k,1}), '--s', 'Color',cols(t,:), 'LineWidth',1.6, 'MarkerSize',4, 'DisplayName', sprintf('%s t-SNE + 6 neighbours', strrep(tags{t},'_',' ')));
            r = q('tsne',0);         plot(ax, r.sweeps, r.(met{k,1}), ':o', 'Color',cols(t,:), 'LineWidth',1.4, 'MarkerSize',4, 'DisplayName', sprintf('%s t-SNE', strrep(tags{t},'_',' ')));
        end
        xlim(ax, [0 40]); xlabel(ax,'number of sweeps'); ylabel(ax, met{k,2});
        if k == size(met,1), legend(ax,'Location','eastoutside','FontSize',7); end
    end
    title(tl, 'Full-NLOS users: hierarchical search vs t-SNE for growing surfaces');
    exportgraphics(f, fullfile(outDir,'sizes_comparison.png'), 'Resolution',110);  close(f);
end

function write_readme(outDir, T, V, val, fit, charts, nU, nz, Ks, Rs, runtime, pop, nNone, train)
    if strcmp(pop,'nolos')
        ttl = sprintf('# Banana experiment: charts from all users WITHOUT LoS, beam selection vs number of sweeps (noise `%s`)', nz.name);
        who = sprintf(['- **Users (banana):** every user without LoS (%d): reflection only (%d, BS link only via a wall: direct-channel power more than 36 dB below free space) ' ...
                       'and full NLOS (%d, no BS-UE channel). Their BS channel is part of the t-SNE input and of the served SNR. The full-NLOS results compare directly with results_aoa/serving_nlos_size (charts from full-NLOS users only) and results_aoa/serving_sweeps (charts from all users).'], nU, nU - nNone, nNone);
    else
        ttl = sprintf('# Full-NLOS users: EMS beam selection vs number of sweeps, surfaces 60..240 (noise `%s`)', nz.name);
        who = sprintf('- **Users:** full NLOS = no BS-UE channel at all (%d users); the EMSs are their only link.', nU);
    end
    L = {ttl, '', ...
         sprintf('Generated %s by `experiments/Main_AoA_Serving_NLoS_Size.m` (branch aoa-serving), runtime %.0f min.', ...
                 char(datetime('now','Format','yyyy-MM-dd HH:mm')), runtime), '', ...
         '## Setup', ...
         who, ...
         sprintf('- **t-SNE:** one 1-D chart per EMS with the users of this population it reaches (EMS1: %d users, %d anchors, perplexity %d; EMS2: %d, %d, %d). Scored: reachable by at least one EMS and anchor in neither chart.', ...
                 numel(charts(1).users), charts(1).nA, charts(1).perp, numel(charts(2).users), charts(2).nA, charts(2).perp), ...
         '- **Surfaces:** M x M elements, lambda/2, steering along one axis only (vertical: same phase on every column). `M60_raytraced` = the dataset; `M60/120/180/240` = synthesized with spherical waves from the UE, EMS and BS positions, per-user strengths fitted on the 60 x 60 data. Codebook from the beamwidth: N = 68/136/204/272.', ...
         sprintf('- **Schemes vs number of sweeps:** sequential, spread, hierarchical (4 spread entries, then halving), spread + R neighbours, t-SNE AoA + LUT (K = %s spread training entries; %s), t-SNE + LUT entry + R neighbours (R = %s), exhaustive (N), true AoA + LUT.', mat2str(Ks), trtxt(train), mat2str(Rs)), ...
         '- Measurement-based selection uses each EMS path''s own gain (noise-free, in favour of the sweeping baselines); the served SNR uses the true channel; 23 dBm, -92 dBm.', '', ...
         '## Validation of the synthesized channels (60 x 60 vs ray-traced)', ''};
    for i = 1:size(val,1)
        L{end+1} = sprintf('- EMS%d: element grid orientation (m, n) = (%+d, %+d), phase sign %+d; |correlation| of the EMS->BS channel with the model %.3f; UE->EMS channels: median |correlation| %.3f, above 0.9 for %.0f%% of the %d reachable full-NLOS users.', ...
            val{i,1}, fit(val{i,1}).ori(1), fit(val{i,1}).ori(2), fit(val{i,1}).ori(3), val{i,2}, val{i,3}, val{i,4}, val{i,5}); %#ok<AGROW>
    end
    vv = V.note(startsWith(V.note,'validation'));
    for i = 1:numel(vv), L{end+1} = ['- ' vv{i}]; end %#ok<AGROW>
    L = [L, {'', '## Median SNR gain vs zero phase [dB] / 10th percentile [dB], full-NLOS users', ''}];
    S = [1 2 3 4 6 8 10 12 16 20];
    for sz = unique(T.surface,'stable').'
        Tg = T(strcmp(T.surface,sz{1}) & strcmp(T.group,'full NLOS, reachable'),:);
        L{end+1} = sprintf('**%s** (codebook %d, %d users)', sz{1}, Tg.N(1), Tg.n_users(1)); %#ok<AGROW>
        L{end+1} = ['| scheme | ' strjoin(arrayfun(@num2str, S, 'UniformOutput', false), ' | ') ' |']; %#ok<AGROW>
        L{end+1} = ['|---|' repmat('---|', 1, numel(S))]; %#ok<AGROW>
        for sc = {'sequential',0; 'spread',0; 'hierarchical',0; 'spread+refine',6; 'tsne',0; 'tsne+refine',2; 'tsne+refine',6}.'
            c = arrayfun(@(s) cellv(Tg, sc{1}, sc{2}, s), S, 'UniformOutput', false);
            nm = sc{1};  if sc{2} > 0, nm = sprintf('%s (R=%d)', nm, sc{2}); end
            L{end+1} = sprintf('| %s | %s |', nm, strjoin(c, ' | ')); %#ok<AGROW>
        end
        for sc = {'true AoA + LUT','exhaustive'}
            r = Tg(strcmp(Tg.scheme, sc{1}),:);
            L{end+1} = sprintf('| %s (%d sweeps) | %.1f / %.1f |%s', sc{1}, r.sweeps, r.dsnr_median_dB, r.dsnr_p10_dB, repmat(' |', 1, numel(S)-1)); %#ok<AGROW>
        end
        r = Tg(strcmp(Tg.scheme,'baseline'),:);
        L{end+1} = sprintf('Zero-phase median SNR %.1f dB, exhaustive median SNR %.1f dB.', r.snr_median_dB, Tg.snr_median_dB(strcmp(Tg.scheme,'exhaustive'))); %#ok<AGROW>
        vv = V.note(strcmp(V.surface, sz{1}) & startsWith(V.note,'K='));
        L{end+1} = ['t-SNE angle error / LUT: ' strjoin(vv, '; ')]; %#ok<AGROW>
        L{end+1} = ''; %#ok<AGROW>
    end
    L{end+1} = 'Figures: `sweeps_<surface>.png` per surface, `sizes_comparison.png` across surfaces. All numbers: `summary.csv`, `validation_and_aoa.csv`.';
    fid = fopen(fullfile(outDir,'README.md'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
end

function c = cellv(T, sc, R, s)
    r = T(strcmp(T.scheme,sc) & T.R==R & T.sweeps==s,:);
    if isempty(r), c = '-'; else, c = sprintf('%.1f / %.1f', r.dsnr_median_dB(1), r.dsnr_p10_dB(1)); end
end

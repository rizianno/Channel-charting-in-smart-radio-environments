%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EMS BEAM SELECTION vs NUMBER OF SWEEPS: t-SNE AoA (+ refinement) vs beam sweeping
%
% Every scheme chooses one serving codebook entry per EMS for each user; the result is
% shown against the number of sweeps (acquisitions) it needed, with the codebook fixed.
%   sequential      : entries 1, 2, 3, ... in order, best so far            (worst baseline)
%   spread          : s entries spread evenly over the codebook, best of them
%   hierarchical    : 4 spread entries, then the two neighbours at half the spacing around the
%                     current best, repeated down to spacing 1 (coarse-to-fine)
%   tsne            : K spread training sweeps -> t-SNE azimuth -> look-up table      (K sweeps)
%   tsne+R          : tsne, then measure the LUT entry and R neighbours (R/2 per side),
%                     keep the best                                            (K + R + 1 sweeps)
%   spread+R        : B spread sweeps, then R neighbours (R/2 per side) of the best (B + R)
%   exhaustive      : all N entries per EMS (N sweeps), and all N^2 pairs (joint upper bound)
% Selection by sweeping uses the measured gain of each EMS path separately: each EMS reaches
% the BS from one fixed direction (its EMS->BS channel is rank one), so the BS separates the two
% surfaces spatially and both are swept in parallel; measurements are taken noise-free (in
% favour of the sweeping baselines). The served SNR always uses the full true channel.
%
% Codebook: N entries, phase step (i - N/2) * 360/N deg/element over -180..180. N follows from
% the beamwidth: the half-power width of the steered beam in sin(theta) is 0.886 lambda/D with
% D = M lambda/2 (lambda/2 spacing), and N = ceil(2 / width) covers sin(theta) in [-1, 1] with
% beams crossing at -3 dB (M = 60 -> N = 68).
% Larger surface: EMS gains scaled by SRV_GAIN_DB (19.1 dB = 180x180 vs 60x60 elements,
% coherent gain ~ (number of elements)^2); the beams stay those of the 60x60 data (a true
% larger aperture needs new channels).
% t-SNE AoA: joint chart, all users in t-SNE, 15% of the users in front of both surfaces as
% anchors, one chart per training sweep (both EMSs on the same entry), label-free local fusion;
% LUT: steering law sin(theta) = sin(theta_specular) + step/180 with the known BS direction.
% Link: uplink, 23 dBm, -92 dBm, SNR = P/N0 * lambda_max(H'H). Users: unlabelled, grouped by
% the BS-UE link (LoS / reflection only / none, LoS from the direct-channel power vs free
% space, threshold -36 dB).
%
% Environment: SRV_NOISE (est_G60dB), SRV_GAIN_DB ("0,19.1"), SRV_KS ("1:8"), SRV_RS ("2,4,6"),
%              SRV_SMAX (40), SRV_WORKERS (48), SRV_MAX_USERS (testing)
% Outputs: results_aoa/serving_sweeps/<noise>/ (README.md, summary.csv, figures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();  params.target = 'aoa';
CAT  = abl_catalog();
nz0  = CAT.noise(strcmp({CAT.noise.name}, envs('SRV_NOISE','est_G60dB')));
gdBs = str2double(strsplit(envs('SRV_GAIN_DB','0,19.1'), ','));
Ks   = eval(envs('SRV_KS','1:8'));
Rs   = str2double(strsplit(envs('SRV_RS','2,4,6'), ','));
Smax = str2double(envs('SRV_SMAX','40'));
nW   = str2double(envs('SRV_WORKERS','48'));

M = params.M;                                   % elements along the steering axis
width_u = 0.886 * 2 / M;                        % half-power beamwidth in sin(theta), lambda/2
N = ceil(2 / width_u);  N = N + mod(N,2);       % even: entry N/2 = zero phase
step = ((1:N) - N/2) * 360/N;
spread = @(s) unique(round((2*(1:s) - 1) * N / (2*s)));
fprintf('Beamwidth (sin) %.4f -> codebook N = %d entries, %.2f deg/element\n', width_u, N, 360/N);

%% ------------------------------------------------------------------------
% DATA, USERS, LINK CLASSES
%% ------------------------------------------------------------------------
Data = load_dataset(params);
maxU = str2double(envs('SRV_MAX_USERS',''));
if ~isnan(maxU)
    for fld = {'R_xy','H_LOS_1bounce','H_I_RIS1','H_I_RIS2'}
        x = Data.(fld{1});  Data.(fld{1}) = x(1:maxU,:,:);
    end
end
U  = abl_users(Data, params);
ch = abl_charts(U, 'S3');
lab = false(U.n,1);  lab(ch.users(1:ch.nA)) = true;
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);
lam = 3e8/30e9;
Pd  = squeeze(sum(sum(abs(Data.H_LOS_1bounce).^2,2),3)) / (Nb*Nr);
dd  = sqrt(sum(([U.xy, params.UE_height*ones(U.n,1)] - U.BS_position).^2, 2));
los  = U.has_direct & 10*log10(Pd ./ (lam./(4*pi*dd)).^2) > -36;
refl = U.has_direct & ~los;
none = ~U.has_direct;
rho  = 10^((nz0.tx_dBm - nz0.noise_dBm)/10);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);
specular = -azBS;
Hd = reshape(permute(double(Data.H_LOS_1bounce),[2 3 1]), Nb*Nr, U.n);
lmax = @(H) lam_max(H, Nb, Nr);

pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', nW); end
if ~exist(fullfile(root,'cache_ablations'),'dir'), mkdir(fullfile(root,'cache_ablations')); end
rows = {};  curves = struct();
for gdB = gdBs
    gtag = sprintf('ems%+.1fdB', gdB);
    G = 10^(gdB/20);                                  % amplitude scaling of the EMS paths
    Dg = Data;  Dg.H_I_RIS1 = Data.H_I_RIS1 * G;  Dg.H_I_RIS2 = Data.H_I_RIS2 * G;
    nz = nz0;  if gdB ~= 0, nz.name = sprintf('%s_%s', nz0.name, gtag); end

    % EMS paths of every entry, and their own (separately measured) gains
    A = cell(1,2);  Gown = cell(1,2);
    for e = 1:2
        Ho = double(Dg.(sprintf('H_o_RIS%d',e)));
        X  = reshape(permute(Dg.(sprintf('H_I_RIS%d',e)),[2 3 1]), params.M*params.N, []);
        A{e} = zeros(Nb*Nr, U.n, N);  Gown{e} = zeros(U.n, N);
        for i = 1:N
            th = exp(-1j*kron(deg2rad(step(i))*(1:params.M), ones(1,params.N))).';
            A{e}(:,:,i) = reshape(Ho*(th.*X), Nb*Nr, U.n);
            Gown{e}(:,i) = lmax(A{e}(:,:,i)).';
        end
    end
    snrOf = @(i1, i2) rho * lmax(Hd + pick(A{1}, i1) + pick(A{2}, i2));
    base  = snrOf(repmat(N/2,U.n,1), repmat(N/2,U.n,1));
    jbest = zeros(1,U.n);
    for i1 = 1:N, for i2 = 1:N, jbest = max(jbest, lmax(Hd + A{1}(:,:,i1) + A{2}(:,:,i2))); end, end
    jointEx = rho * jbest;
    ok  = ~lab & (jointEx > 0).';                     % unlabelled users with any channel
    grp = {'no LoS (reflection only or no channel)', ok & ~los; 'reflection only', ok & refl; ...
           'no BS-UE channel', ok & none; 'LoS', ok & los; 'all users', ok};
    sel = @(cand) best_of(Gown, cand);               % cand{e}: U x c entry indices

    add = @(rows, scheme, R, s, snr) [rows; collect(gdB, scheme, R, s, snr, base, grp)];
    rows = add(rows, 'exhaustive_joint', 0, N^2, jointEx);
    [i1, i2] = sel({repmat(1:N,U.n,1), repmat(1:N,U.n,1)});
    rows = add(rows, 'exhaustive', 0, N, snrOf(i1, i2));
    rows = add(rows, 'baseline', 0, 0, base);

    % ---- sweeping baselines, anytime curves
    for s = 1:Smax
        [i1, i2] = sel({repmat(1:s,U.n,1), repmat(1:s,U.n,1)});
        rows = add(rows, 'sequential', 0, s, snrOf(i1, i2));
        sp = spread(s);
        [i1, i2] = sel({repmat(sp,U.n,1), repmat(sp,U.n,1)});
        rows = add(rows, 'spread', 0, s, snrOf(i1, i2));
        [i1, i2] = hierarchical(Gown, N, s, spread);
        rows = add(rows, 'hierarchical', 0, s, snrOf(i1, i2));
        for R = Rs
            if s - R < 1, continue; end
            sp = spread(s - R);
            [b1, b2] = sel({repmat(sp,U.n,1), repmat(sp,U.n,1)});
            [i1, i2] = sel({nbr(b1, R, N), nbr(b2, R, N)});
            rows = add(rows, 'spread+refine', R, s, snrOf(i1, i2));
        end
    end

    % ---- t-SNE training charts (cached)
    cacheFile = fullfile(root,'cache_ablations', sprintf('serving_%s_U%d.mat', nz.name, U.n));
    Rmap = load_cache(cacheFile);
    P = unique(cell2mat(arrayfun(@(K) step(spread(K)).', Ks, 'UniformOutput', false).'), 'stable');
    P = [P P];
    keys_ = arrayfun(@(i) sprintf('%s|%s|F0|%s', nz.name, ch.id, abl_stkey(true, P(i,:))), 1:size(P,1), 'UniformOutput', false);
    todo = find(~cellfun(@(k) isKey(Rmap,k), keys_));
    fprintf('\n[%s] %d training sweeps, %d charts to compute\n', nz.name, size(P,1), numel(todo));
    futs = parallel.FevalFuture.empty;  fk = {};  nx = 1;  done = 0;
    while nx <= numel(todo) || ~isempty(futs)
        while nx <= numel(todo) && numel(futs) < 2*nW
            i = todo(nx);  nx = nx + 1;
            st = struct('direct', true, 'phi', P(i,:), 'key', abl_stkey(true, P(i,:)));
            F = abl_state_features(Dg, params, st, nz);
            D = squareform(pdist(F(ch.users,:)));
            futs(end+1) = parfeval(pool, @abl_chart, 1, D, ch.labels, ch.nA, ch.valid, ...
                                   params.tsne_perplexity, params.K_fraction); %#ok<SAGROW>
            fk{end+1} = keys_{i}; %#ok<SAGROW>
        end
        [j, out] = fetchNext(futs);
        Rmap(fk{j}) = out;  futs(j) = [];  fk(j) = [];  done = done + 1;
        if mod(done,10) == 0 || (nx > numel(todo) && isempty(futs))
            save_cache(cacheFile, Rmap);
            fprintf('  %d/%d charts (%.1f min)\n', done, numel(todo), toc(t0)/60);
        end
    end

    % ---- t-SNE AoA -> LUT (-> refinement)
    lutT = lut_physics(U.az, specular, step);
    rows = add(rows, 'true AoA + LUT', 0, 0, snrOf(lutT(:,1), lutT(:,2)));
    for K = Ks
        Pk = step(spread(K)).';  Pk = [Pk Pk];
        outs = Rmap(sprintf('%s|%s|F0|%s', nz.name, ch.id, abl_stkey(true, Pk(1,:))));
        for s = 2:size(Pk,1), outs(s) = Rmap(sprintf('%s|%s|F0|%s', nz.name, ch.id, abl_stkey(true, Pk(s,:)))); end
        E3 = cat(3, outs.est);
        if size(Pk,1) == 1, est = E3; else, [~, est] = abl_fuse(E3, [outs.q_free], ch.nA, ch.labels); end
        azHat = nan(U.n,2);  azHat(ch.users,:) = est;
        lut = lut_physics(azHat, specular, step);
        nK  = size(Pk,1);
        rows = add(rows, 'tsne', 0, nK, snrOf(lut(:,1), lut(:,2)));
        for R = Rs
            [i1, i2] = sel({nbr(lut(:,1), R, N, true), nbr(lut(:,2), R, N, true)});
            rows = add(rows, 'tsne+refine', R, nK + R + 1, snrOf(i1, i2));
        end
        err = mod(azHat - U.az + 180, 360) - 180;  fb = U.both & ~lab;
        fprintf('  K=%d: AoA RMSE %.1f/%.1f deg; LUT entry off by >1 from the true-AoA LUT: %.0f%%/%.0f%%\n', ...
                nK, sqrt(mean(err(fb,1).^2)), sqrt(mean(err(fb,2).^2)), ...
                100*mean(cdist(lut(fb,1),lutT(fb,1),N) > 1), 100*mean(cdist(lut(fb,2),lutT(fb,2),N) > 1));
    end
end

T = cell2table(rows, 'VariableNames', {'gain_dB','scheme','R','sweeps','group','n_users', ...
    'dsnr_median_dB','dsnr_mean_dB','snr_median_dB','se_mean','dse_mean','coverage_0dB','coverage_m10dB'});
outDir = fullfile(root,'results_aoa','serving_sweeps', nz0.name);
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir,'summary.csv'));
for gdB = gdBs, plot_sweeps(T(T.gain_dB == gdB,:), gdB, N, Rs, outDir); end
write_readme(outDir, T, gdBs, N, step, width_u, M, Ks, Rs, Smax, nz0, U, ch, los, refl, none, specular, toc(t0)/60);
fprintf('\n=== DONE (%.1f min) ===\n', toc(t0)/60);


%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function l = lam_max(H, Nb, Nr)
    H  = reshape(H, Nb, Nr, []);
    h1 = squeeze(H(:,1,:));  h2 = squeeze(H(:,2,:));
    g11 = sum(abs(h1).^2,1);  g22 = sum(abs(h2).^2,1);  g12 = sum(conj(h1).*h2,1);
    l = (g11+g22)/2 + sqrt(((g11-g22)/2).^2 + abs(g12).^2);
end

function H = pick(Ae, idx)
    [d, n, ~] = size(Ae);
    H = reshape(Ae(:, (idx(:).'-1)*n + (1:n)), d, n);
end

function [i1, i2] = best_of(Gown, cand)
% per user and EMS: the candidate entry with the largest separately measured EMS-path gain
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
% R neighbours (R/2 on each side, circular) of entry b, optionally plus b itself
    if nargin < 4, withCentre = false; end
    off = [-(R/2):-1, 1:(R/2)];
    if withCentre, off = [0 off]; end
    c = mod(b(:) - 1 + off, N) + 1;
end

function d = cdist(a, b, N)
    d = abs(a - b);  d = min(d, N - d);
end

function [i1, i2] = hierarchical(Gown, N, s, spread)
% coarse-to-fine search per EMS with s sweeps: 4 spread entries, then the two neighbours at
% half the spacing around the current best (one at a time when s is odd)
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

function sel = lut_physics(az, specular, step)
    sel = zeros(size(az,1), 2);
    for e = 1:2
        want = mod(180*(sind(az(:,e)) - sind(specular(e))) + 180, 360) - 180;
        [~, sel(:,e)] = min(abs(mod(want - step + 180, 360) - 180), [], 2);
        sel(isnan(az(:,e)), e) = numel(step)/2;
    end
end

function r = collect(gdB, scheme, R, s, snr, base, grp)
    r = cell(size(grp,1), 13);
    for g = 1:size(grp,1)
        m = grp{g,2};  x = snr(m);  b = base(m);
        d = 10*log10(x ./ b);  d = d(isfinite(d));
        r(g,:) = {gdB, scheme, R, s, grp{g,1}, nnz(m), median(d), mean(d), median(10*log10(x)), ...
                  mean(log2(1+x)), mean(log2(1+x) - log2(1+b)), 100*mean(x >= 1), 100*mean(x >= 0.1)};
    end
end

function R = load_cache(file)
    R = containers.Map();
    if isfile(file)
        C = load(file);
        for i = 1:numel(C.keys), R(C.keys{i}) = C.vals{i}; end
    end
end

function save_cache(file, R)
    keys = R.keys;  vals = R.values; %#ok<NASGU>
    save(file, 'keys', 'vals', '-v7.3');
end

function plot_sweeps(T, gdB, N, Rs, outDir)
    groups = {'no LoS (reflection only or no channel)', 'reflection only', 'no BS-UE channel'};
    met = {'dsnr_median_dB','median SNR gain vs zero phase [dB]'; 'dsnr_mean_dB','mean SNR gain vs zero phase [dB]'; ...
           'dse_mean','mean SE gain [bit/s/Hz]'};
    blues = [0.30 0.75 0.93; 0.00 0.45 0.74; 0.00 0.20 0.45];
    orng  = [1.00 0.75 0.40; 0.95 0.55 0.15; 0.70 0.30 0.00];
    f = figure('Position',[50 50 1800 1100]);
    tl = tiledlayout(f, numel(groups), size(met,1), 'TileSpacing','compact');
    for g = 1:numel(groups)
        for k = 1:size(met,1)
            ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
            q = @(sc, R) sortrows(T(strcmp(T.scheme,sc) & T.R==R & strcmp(T.group,groups{g}),:), 'sweeps');
            r = q('sequential',0);   plot(ax, r.sweeps, r.(met{k,1}), ':', 'Color',[0.5 0.5 0.5], 'LineWidth',1.5, 'DisplayName','sequential sweep');
            r = q('spread',0);       plot(ax, r.sweeps, r.(met{k,1}), '-', 'Color',[0.85 0.33 0.10], 'LineWidth',1.8, 'DisplayName','spread sweep');
            r = q('hierarchical',0); plot(ax, r.sweeps, r.(met{k,1}), '-', 'Color',[0.49 0.18 0.56], 'LineWidth',1.8, 'DisplayName','hierarchical search');
            for i = 1:numel(Rs)
                r = q('spread+refine',Rs(i)); plot(ax, r.sweeps, r.(met{k,1}), '--', 'Color',orng(i,:), 'LineWidth',1.2, ...
                    'DisplayName', sprintf('spread sweep + %d neighbours', Rs(i)));
            end
            r = q('tsne',0); plot(ax, r.sweeps, r.(met{k,1}), '-o', 'Color',[0 0.45 0.74], 'LineWidth',2, 'MarkerSize',5, ...
                    'MarkerFaceColor',[0 0.45 0.74], 'DisplayName','t-SNE AoA + LUT');
            for i = 1:numel(Rs)
                r = q('tsne+refine',Rs(i)); plot(ax, r.sweeps, r.(met{k,1}), '-s', 'Color',blues(i,:), 'LineWidth',1.5, 'MarkerSize',5, ...
                    'DisplayName', sprintf('t-SNE AoA + LUT + %d neighbours', Rs(i)));
            end
            r = q('true AoA + LUT',0);   yline(ax, r.(met{k,1}), ':', 'Color',[0 0.45 0.74], 'LineWidth',1.5, 'DisplayName','true AoA + LUT (no sweeps)');
            r = q('exhaustive',0);       yline(ax, r.(met{k,1}), '--k', 'LineWidth',1.5, 'DisplayName', sprintf('exhaustive (%d sweeps)', N));
            r = q('exhaustive_joint',0); yline(ax, r.(met{k,1}), ':k', 'LineWidth',1.5, 'DisplayName', sprintf('joint exhaustive (%d sweeps)', N^2));
            xlim(ax, [0 max(T.sweeps(T.sweeps <= 60))]);
            xlabel(ax, 'number of sweeps'); ylabel(ax, met{k,2}); title(ax, groups{g});
            if g == 1 && k == size(met,1), legend(ax, 'Location','eastoutside', 'FontSize',8); end
        end
    end
    title(tl, sprintf('EMS beam selection vs number of sweeps, codebook N = %d per EMS, EMS gain %+.1f dB', N, gdB));
    exportgraphics(f, fullfile(outDir, sprintf('sweeps_gain%+.1fdB.png', gdB)), 'Resolution',110);
    close(f);
end

function write_readme(outDir, T, gdBs, N, step, width_u, M, Ks, Rs, Smax, nz, U, ch, los, refl, none, specular, runtime)
    L = {sprintf('# EMS beam selection vs number of sweeps (noise `%s`)', nz.name), '', ...
         sprintf('Generated %s by `experiments/Main_AoA_Serving_Sweeps.m` (branch aoa-serving), runtime %.0f min.', ...
                 char(datetime('now','Format','yyyy-MM-dd HH:mm')), runtime), '', ...
         '## Codebook from the beamwidth', ...
         sprintf('- %d elements along the steering axis at lambda/2 -> aperture %d lambda (%.2f m at 30 GHz). Half-power beamwidth in sin(theta): 0.886 lambda / D = %.4f (%.1f deg near the normal, %.1f deg at 60 deg).', ...
                 M, M/2, M/2*0.01, width_u, asind(width_u), rad2deg(width_u/cosd(60))), ...
         sprintf('- Covering sin(theta) in [-1, 1] with beams crossing at -3 dB: N = %d entries per EMS, phase step %.2f deg/element, over -180..180.', N, 360/N), ...
         sprintf('- 180 elements per side would give width %.4f and about %d entries; here only its gain is emulated (+19.1 dB = (180/60)^4), the beams stay those of the 60-element data.', ...
                 0.886*2/180, ceil(2/(0.886*2/180))), ...
         sprintf('- Far-field distance 2 D^2 / lambda = %.0f m: users closer than that are in the near field of the surface, where the steered (far-field) beams lose gain.', 2*(M/2*0.01)^2/0.01), ...
         sprintf('- Specular direction of EMS1 / EMS2 (BS mirrored): %.1f / %.1f deg.', specular(1), specular(2)), '', ...
         '## Schemes (all vs number of sweeps, codebook fixed)', ...
         '- `sequential`: entries 1, 2, 3, ... (worst baseline). `spread`: s entries spread evenly. `hierarchical`: 4 spread entries, then the two neighbours at half the spacing around the current best, down to spacing 1.', ...
         sprintf('- `tsne`: K spread training sweeps (K = %s, both EMSs on the same entry) -> t-SNE azimuth (joint chart, all %d users, %d anchors, label-free local fusion) -> look-up table (steering law with the known BS direction). `tsne+refine` (R = %s): then measure the LUT entry and R neighbours (R/2 per side): K + R + 1 sweeps. `spread+refine`: B spread sweeps + R neighbours of the best: B + R sweeps.', ...
                 mat2str(Ks), U.n, ch.nA, mat2str(Rs)), ...
         sprintf('- Upper bounds: `exhaustive` (all %d entries per EMS) and `exhaustive_joint` (all %d pairs). `true AoA + LUT`: the look-up table with the exact angle (no sweeps).', N, N^2), ...
         '- Selection by measurement uses each EMS path''s own gain (the BS separates the two surfaces spatially: each EMS->BS channel is rank one), noise-free (in favour of the sweeping baselines); t-SNE uses estimated channels. The served SNR uses the full true channel (direct + both EMSs).', ...
         sprintf('- Users: unlabelled with any channel; LoS %d, reflection only %d, no BS-UE channel %d (all users). Link: 23 dBm, -92 dBm, SNR = P/N0 lambda_max(H''H).', nnz(los), nnz(refl), nnz(none)), ''};
    for gdB = gdBs
        L{end+1} = sprintf('## EMS gain %+.1f dB: median SNR gain vs zero phase [dB], users without LoS', gdB); %#ok<AGROW>
        L{end+1} = ['| sweeps | ' strjoin(arrayfun(@num2str, [1 2 3 4 5 6 8 10 12 16 20 30], 'UniformOutput', false), ' | ') ' |']; %#ok<AGROW>
        L{end+1} = ['|---|' repmat('---|', 1, 12)]; %#ok<AGROW>
        Tg = T(T.gain_dB == gdB & strcmp(T.group,'no LoS (reflection only or no channel)'),:);
        for sc = {'sequential',0; 'spread',0; 'hierarchical',0; 'spread+refine',2; 'spread+refine',4; 'spread+refine',6; ...
                  'tsne',0; 'tsne+refine',2; 'tsne+refine',4; 'tsne+refine',6}.'
            c = arrayfun(@(s) val(Tg, sc{1}, sc{2}, s), [1 2 3 4 5 6 8 10 12 16 20 30], 'UniformOutput', false);
            nm = sc{1};  if sc{2} > 0, nm = sprintf('%s (R=%d)', nm, sc{2}); end
            L{end+1} = sprintf('| %s | %s |', nm, strjoin(c, ' | ')); %#ok<AGROW>
        end
        for sc = {'true AoA + LUT','exhaustive','exhaustive_joint'}
            r = Tg(strcmp(Tg.scheme, sc{1}),:);
            L{end+1} = sprintf('| %s (%d sweeps) | %.1f |%s', sc{1}, r.sweeps, r.dsnr_median_dB, repmat(' |', 1, 11)); %#ok<AGROW>
        end
        L{end+1} = ''; %#ok<AGROW>
    end
    L{end+1} = 'Figures: `sweeps_gain<g>dB.png` (median / mean SNR gain and SE gain vs sweeps for users without LoS, reflection-only and no-channel users). All numbers, also LoS and all users and coverage: `summary.csv`.';
    fid = fopen(fullfile(outDir,'README.md'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
end

function c = val(T, sc, R, s)
    r = T(strcmp(T.scheme,sc) & T.R==R & T.sweeps==s,:);
    if isempty(r), c = '-'; else, c = sprintf('%.1f', r.dsnr_median_dB(1)); end
end

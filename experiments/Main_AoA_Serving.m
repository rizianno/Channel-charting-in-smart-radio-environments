%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SERVING USERS FROM AoA ESTIMATED WITH FEW EMS ACQUISITIONS
%
% End-to-end: K training acquisitions -> semi-supervised t-SNE estimates the azimuth of
% each user seen from EMS1 and EMS2 -> a look-up table maps each estimated angle to one
% entry of a serving codebook -> uplink SNR / spectral efficiency / coverage of the user
% with that EMS configuration, compared with:
%   baseline        : both EMSs with zero phase (no configuration)
%   best_training   : best of the K training configurations measured for that user
%                     (few sweeps, no angle estimation)
%   genie_lut       : look-up table fed with the TRUE angle (cost of the angle error)
%   full_sweep      : best of all N x N serving configurations (exhaustive, upper bound)
%
% Serving codebook: N entries per EMS, phase step between adjacent elements
%   step(i) = (i - N/2) * 360/N deg/element, i = 1..N   (entry N/2 = 0 = baseline)
% For lambda/2 spacing the beam goes to sin(theta) = sin(theta_specular) + step/180.
% Training entries: K entries spread evenly over the codebook, round((2k-1) N / (2K)),
%   e.g. N = 30: K=1 -> 15, K=3 -> 5,15,25 (never neighbouring entries).
%   mode "diag": acquisition k sets both EMSs to training entry k  -> K acquisitions
%   mode "grid": all pairs of training entries                    -> K^2 acquisitions
% AoA estimation (best setup of the ablations): one joint chart (azimuth from EMS1, EMS2),
%   all users in t-SNE, 15% of the users in front of both surfaces as anchors, one chart
%   per acquisition (log-covariance, Log-Euclidean distance, tsne_sem) and local fusion
%   with the LABEL-FREE quality (no ground truth of the served users is used).
% Look-up tables (angle -> entry, per EMS):
%   physics : step = 180 (sin(theta_hat) - sin(theta_specular)), nearest entry; the specular
%             direction follows from the known BS and EMS positions
%   data    : entry that served the 7 anchors with the closest true angle best (majority),
%             i.e. an offline calibration with the labelled users
% Link: uplink, UE power 23 dBm, noise -92 dBm, best rank-1 beamforming (SNR =
%   P/N0 * lambda_max(H'H)), H = direct + EMS1 + EMS2 (true channel). SE = log2(1+SNR).
% Scored users: unlabelled users (anchors excluded), groups all / in front of both
%   surfaces / without a direct path.
%
% Environment: SRV_NOISE ("est_G60dB,noiseless_1e6"), SRV_MODES ("diag,grid"),
%   SRV_KS ("1,2,3,4,7"), SRV_N (30), SRV_WORKERS (48), SRV_MAX_USERS (testing)
% Charts are cached in cache_ablations/serving_<noise>_U<n>.mat (git-ignored); charts
% already computed by Main_AoA_Ablations for the same acquisitions are reused (read-only).
% Outputs: results_aoa/serving/<noise>/ (README.md, summary.csv, figures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();  params.target = 'aoa';
CAT    = abl_catalog();
noises = strtrim(strsplit(envs('SRV_NOISE','est_G60dB,noiseless_1e6'), ','));
modes  = strtrim(strsplit(envs('SRV_MODES','diag,grid'), ','));
Ks     = str2double(strsplit(envs('SRV_KS','1,2,3,4,7'), ','));
N      = str2double(envs('SRV_N','30'));
nW     = str2double(envs('SRV_WORKERS','48'));
assert(mod(N,2) == 0, 'SRV_N must be even (entry N/2 = zero phase)');
step  = ((1:N) - N/2) * 360/N;                     % deg/element
train = @(K) round((2*(1:K) - 1) * N / (2*K));     % training entry indices

%% ------------------------------------------------------------------------
% DATA, USERS, CHART
%% ------------------------------------------------------------------------
Data = load_dataset(params);
maxU = str2double(envs('SRV_MAX_USERS',''));
if ~isnan(maxU)
    for fld = {'R_xy','H_LOS_1bounce','H_I_RIS1','H_I_RIS2'}
        x = Data.(fld{1});  Data.(fld{1}) = x(1:maxU,:,:);
    end
end
U  = abl_users(Data, params);
ch = abl_charts(U, 'S3');                          % joint chart, all users
lab = false(U.n,1);  lab(ch.users(1:ch.nA)) = true;
unl = ~lab;
fprintf('Users %d, anchors %d, codebook N = %d\n', U.n, ch.nA, N);

% EMS-reflected channel of every serving entry, and the direct channel (Nb*Nr x U)
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);
Hd = reshape(permute(double(Data.H_LOS_1bounce),[2 3 1]), Nb*Nr, U.n);
A  = cell(1,2);
for e = 1:2
    Ho = double(Data.(sprintf('H_o_RIS%d',e)));
    Hi = Data.(sprintf('H_I_RIS%d',e));
    X  = reshape(permute(Hi,[2 3 1]), size(Hi,2), []);
    A{e} = zeros(Nb*Nr, U.n, N);
    for i = 1:N
        th = exp(-1j*kron(deg2rad(step(i))*(1:params.M), ones(1,params.N))).';
        A{e}(:,:,i) = reshape(Ho*(th.*X), Nb*Nr, U.n);
    end
end
lmax = @(H) lam_max(H, Nb, Nr);

% exhaustive full sweep (upper bound) and baseline
nz0 = CAT.noise(strcmp({CAT.noise.name}, noises{1}));
rho = 10^((nz0.tx_dBm - nz0.noise_dBm)/10);
best_full = zeros(1,U.n);  arg_full = zeros(2,U.n);
for i1 = 1:N
    for i2 = 1:N
        l = lmax(Hd + A{1}(:,:,i1) + A{2}(:,:,i2));
        b = l > best_full;  best_full(b) = l(b);  arg_full(:,b) = repmat([i1;i2],1,nnz(b));
    end
end
snr_base = rho * lmax(Hd + A{1}(:,:,N/2) + A{2}(:,:,N/2));
snr_full = rho * best_full;

% EMS-path gain of every entry (for the data look-up table and its calibration)
gainE = cellfun(@(Ae) squeeze(lam_max_stack(Ae, Nb, Nr)), A, 'UniformOutput', false);   % U x N
[~, bestEntry] = cellfun(@(g) max(g, [], 2), gainE, 'UniformOutput', false);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);   % BS seen from each EMS
specular = -azBS;                                                          % mirror direction

%% ------------------------------------------------------------------------
% CHARTS OF THE TRAINING ACQUISITIONS (parallel, cached)
%% ------------------------------------------------------------------------
pool = gcp('nocreate');  if isempty(pool), pool = parpool('Processes', nW); end
if ~exist(fullfile(root,'cache_ablations'),'dir'), mkdir(fullfile(root,'cache_ablations')); end
for nzName = noises
    nz = CAT.noise(strcmp({CAT.noise.name}, nzName{1}));
    cacheFile = fullfile(root,'cache_ablations', sprintf('serving_%s_U%d.mat', nz.name, U.n));
    R = load_cache({cacheFile, fullfile(root,'cache_ablations', sprintf('%s_U%d.mat', nz.name, U.n))});
    P = zeros(0,2);
    for m = modes, for K = Ks, P = [P; acq(m{1}, step(train(K)))]; end, end %#ok<AGROW>
    P = unique(P, 'rows', 'stable');
    keys_ = arrayfun(@(i) sprintf('%s|%s|F0|%s', nz.name, ch.id, abl_stkey(true, P(i,:))), 1:size(P,1), 'UniformOutput', false);
    todo = find(~cellfun(@(k) isKey(R,k), keys_));
    fprintf('\n[%s] %d training acquisitions, %d charts to compute\n', nz.name, size(P,1), numel(todo));
    Fmap = containers.Map();
    futs = parallel.FevalFuture.empty;  fk = {};  nx = 1;  done = 0;
    while nx <= numel(todo) || ~isempty(futs)
        while nx <= numel(todo) && numel(futs) < 2*nW
            i = todo(nx);  nx = nx + 1;
            st = struct('direct', true, 'phi', P(i,:), 'key', abl_stkey(true, P(i,:)));
            if ~isKey(Fmap, st.key), Fmap(st.key) = abl_state_features(Data, params, st, nz); end
            F = Fmap(st.key);  remove(Fmap, st.key);
            D = squareform(pdist(F(ch.users,:)));
            futs(end+1) = parfeval(pool, @abl_chart, 1, D, ch.labels, ch.nA, ch.valid, ...
                                   params.tsne_perplexity, params.K_fraction); %#ok<SAGROW>
            fk{end+1} = keys_{i}; %#ok<SAGROW>
        end
        [j, out] = fetchNext(futs);
        R(fk{j}) = out;  futs(j) = [];  fk(j) = [];  done = done + 1;
        if mod(done,10) == 0 || (nx > numel(todo) && isempty(futs))
            save_cache(cacheFile, R);
            fprintf('  %d/%d charts (%.1f min)\n', done, numel(todo), toc(t0)/60);
        end
    end

    %% --------------------------------------------------------------------
    % SERVING PER TRAINING MODE AND K
    %% --------------------------------------------------------------------
    rows = {};  cdfData = struct();
    for m = modes
        for K = Ks
            Pk = acq(m{1}, step(train(K)));
            outs = R(sprintf('%s|%s|F0|%s', nz.name, ch.id, abl_stkey(true, Pk(1,:))));
            for s = 2:size(Pk,1), outs(s) = R(sprintf('%s|%s|F0|%s', nz.name, ch.id, abl_stkey(true, Pk(s,:)))); end
            E3 = cat(3, outs.est);
            if size(Pk,1) == 1, est = E3; else, [~, est] = abl_fuse(E3, [outs.q_free], ch.nA, ch.labels); end
            azHat = nan(U.n,2);  azHat(ch.users,:) = est;

            sel.proposed_physics = lut_physics(azHat, specular, step);
            sel.proposed_data    = lut_data(azHat, U, lab, bestEntry, N);
            sel.genie_physics    = lut_physics(U.az, specular, step);
            sel.genie_data       = lut_data(U.az, U, lab, bestEntry, N);
            snr = struct('baseline', snr_base, 'full_sweep', snr_full);
            for f = fieldnames(sel).'
                snr.(f{1}) = rho * lmax(Hd + pick(A{1}, sel.(f{1})(:,1)) + pick(A{2}, sel.(f{1})(:,2)));
            end
            bt = zeros(1,U.n);
            for s = 1:size(Pk,1)
                [~, i1] = min(abs(step - Pk(s,1)));  [~, i2] = min(abs(step - Pk(s,2)));
                bt = max(bt, lmax(Hd + A{1}(:,:,i1) + A{2}(:,:,i2)));
            end
            snr.best_training = rho * bt;

            err = mod(azHat - U.az + 180, 360) - 180;
            fb  = U.both & unl;
            aoa = [sqrt(mean(err(fb,1).^2)), sqrt(mean(err(fb,2).^2)), median(abs(err(fb,1))), median(abs(err(fb,2)))];
            groups = {'all users', unl; 'in front of both surfaces', U.both & unl; 'without direct path', ~U.has_direct & unl};
            for g = 1:size(groups,1)
                gm = groups{g,2};
                for f = fieldnames(snr).'
                    rows(end+1,:) = [{nz.name, m{1}, K, size(Pk,1), f{1}, groups{g,1}, nnz(gm)}, ...
                                     num2cell(link_stats(snr.(f{1})(gm), snr_base(gm))), num2cell(aoa)]; %#ok<AGROW>
                end
            end
            if K == 3, cdfData.(m{1}) = snr; end
            fprintf('  %s K=%d: median dSNR proposed(physics) %.1f dB, data %.1f dB, best-of-training %.1f dB, full sweep %.1f dB\n', ...
                    m{1}, K, median(10*log10(snr.proposed_physics(unl)./snr_base(unl))), ...
                    median(10*log10(snr.proposed_data(unl)./snr_base(unl))), ...
                    median(10*log10(snr.best_training(unl)./snr_base(unl))), median(10*log10(snr_full(unl)./snr_base(unl))));
        end
    end
    T = cell2table(rows, 'VariableNames', {'noise','mode','K','acquisitions','scheme','group','n_users', ...
        'snr_median_dB','snr_mean_dB','se_mean','coverage_0dB','coverage_10dB', ...
        'dsnr_median_dB','dsnr_mean_dB','dsnr_p10_dB','dse_mean', ...
        'aoa_rmse_ems1','aoa_rmse_ems2','aoa_median_ems1','aoa_median_ems2'});
    outDir = fullfile(root,'results_aoa','serving',nz.name);
    if ~exist(outDir,'dir'), mkdir(outDir); end
    writetable(T, fullfile(outDir,'summary.csv'));
    plot_trends(T, modes, outDir, nz.name);
    plot_cdfs(cdfData, groups, outDir, nz.name);
    write_readme(outDir, nz, T, modes, Ks, N, step, train, U, ch, specular, toc(t0)/60);
end
fprintf('\n=== DONE (%.1f min) ===\n', toc(t0)/60);


%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function P = acq(mode, s)
% EMS configurations [step1 step2] of the training acquisitions
    if strcmp(mode,'diag'), P = [s(:) s(:)];
    else, [a,b] = ndgrid(s, s); P = [a(:) b(:)]; end
end

function l = lam_max(H, Nb, Nr)
% largest eigenvalue of H'H per user, H given as (Nb*Nr) x U (Nr = 2 UE antennas)
    H  = reshape(H, Nb, Nr, []);
    h1 = squeeze(H(:,1,:));  h2 = squeeze(H(:,2,:));
    g11 = sum(abs(h1).^2,1);  g22 = sum(abs(h2).^2,1);  g12 = sum(conj(h1).*h2,1);
    l = (g11+g22)/2 + sqrt(((g11-g22)/2).^2 + abs(g12).^2);
end

function G = lam_max_stack(Ae, Nb, Nr)
% EMS-path gain of every entry: U x N
    G = zeros(size(Ae,2), size(Ae,3));
    for i = 1:size(Ae,3), G(:,i) = lam_max(Ae(:,:,i), Nb, Nr).'; end
end

function H = pick(Ae, idx)
% per-user entry: column u of Ae(:,:,idx(u))
    [d, n, ~] = size(Ae);
    H = reshape(Ae(:, (idx(:).'-1)*n + (1:n)), d, n);
end

function sel = lut_physics(az, specular, step)
% angle -> entry with the reflection-steering law (lambda/2): step = 180 (sin az - sin spec)
    sel = zeros(size(az,1), 2);
    for e = 1:2
        want = mod(180*(sind(az(:,e)) - sind(specular(e))) + 180, 360) - 180;
        d = abs(mod(want - step + 180, 360) - 180);
        [~, sel(:,e)] = min(d, [], 2);
        sel(isnan(az(:,e)), e) = numel(step)/2;
    end
end

function sel = lut_data(az, U, lab, bestEntry, N)
% angle -> majority best entry of the 7 anchors with the closest true angle (per EMS)
    sel = zeros(size(az,1), 2);
    for e = 1:2
        a = find(lab & U.front(:,e));
        ok = ~isnan(az(:,e));
        nn = knnsearch(U.az(a,e), az(ok,e), 'K', 7);
        be = bestEntry{e}(a);
        sel(ok,e)  = mode(be(nn), 2);
        sel(~ok,e) = N/2;
    end
end

function s = link_stats(snr, base)
    d = 10*log10(snr ./ base);
    s = [median(10*log10(snr)), mean(10*log10(snr)), mean(log2(1+snr)), 100*mean(snr >= 1), ...
         100*mean(snr >= 10), median(d), mean(d), prctile(d,10), mean(log2(1+snr) - log2(1+base))];
end

function R = load_cache(files)
    R = containers.Map();
    for f = files
        if isfile(f{1})
            C = load(f{1});
            for i = 1:numel(C.keys), if ~isKey(R, C.keys{i}), R(C.keys{i}) = C.vals{i}; end, end
        end
    end
end

function save_cache(file, R)
    keys = R.keys;  vals = R.values; %#ok<NASGU>
    save(file, 'keys', 'vals', '-v7.3');
end

function plot_trends(T, modes, outDir, nzName)
    sch = {'proposed_physics','t-SNE AoA + physics LUT',[0 0.45 0.74],'-o'; ...
           'proposed_data','t-SNE AoA + data LUT',[0.30 0.75 0.93],'-o'; ...
           'best_training','best of the training acquisitions',[0.85 0.33 0.10],'-s'; ...
           'genie_physics','true AoA + physics LUT',[0 0.45 0.74],':'; ...
           'genie_data','true AoA + data LUT',[0.30 0.75 0.93],':'; ...
           'full_sweep','full sweep (all N^2 configurations)',[0 0 0],'--'; ...
           'baseline','baseline (zero phase)',[0.5 0.5 0.5],'-'};
    met = {'dsnr_median_dB','median SNR gain vs baseline [dB]'; 'dsnr_mean_dB','mean SNR gain vs baseline [dB]'; ...
           'se_mean','mean spectral efficiency [bit/s/Hz]'; 'coverage_0dB','coverage, SNR \geq 0 dB [%]'};
    grp = unique(T.group, 'stable');
    for m = modes
        f = figure('Position',[50 50 1800 950]);
        tl = tiledlayout(f, numel(grp), size(met,1), 'TileSpacing','compact');
        for g = 1:numel(grp)
            for k = 1:size(met,1)
                ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
                for s = 1:size(sch,1)
                    q = sortrows(T(strcmp(T.mode,m{1}) & strcmp(T.group,grp{g}) & strcmp(T.scheme,sch{s,1}),:), 'acquisitions');
                    plot(ax, q.acquisitions, q.(met{k,1}), sch{s,4}, 'Color',sch{s,3}, 'LineWidth',1.5, ...
                         'MarkerSize',4, 'DisplayName', sch{s,2});
                end
                set(ax,'XScale','log'); xticks(ax, unique(T.acquisitions(strcmp(T.mode,m{1}))));
                xlabel(ax,'training acquisitions'); ylabel(ax, met{k,2});
                title(ax, sprintf('%s', grp{g}));
                if g == 1 && k == size(met,1), legend(ax,'Location','eastoutside','FontSize',8); end
            end
        end
        title(tl, sprintf('Serving with few acquisitions (%s training, noise %s)', m{1}, strrep(nzName,'_','\_')));
        exportgraphics(f, fullfile(outDir, sprintf('serving_vs_acquisitions_%s.png', m{1})), 'Resolution',110);
        close(f);
    end
end

function plot_cdfs(C, groups, outDir, nzName)
    sch = {'baseline','baseline (zero phase)',[0.5 0.5 0.5],'-'; ...
           'best_training','best of the 3 training acquisitions',[0.85 0.33 0.10],'-'; ...
           'proposed_physics','t-SNE AoA + physics LUT',[0 0.45 0.74],'-'; ...
           'proposed_data','t-SNE AoA + data LUT',[0.30 0.75 0.93],'-'; ...
           'genie_physics','true AoA + physics LUT',[0 0.45 0.74],':'; ...
           'full_sweep','full sweep',[0 0 0],'--'};
    for m = fieldnames(C).'
        f = figure('Position',[50 50 1500 450]);
        tl = tiledlayout(f, 1, size(groups,1), 'TileSpacing','compact');
        for g = 1:size(groups,1)
            ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
            for s = 1:size(sch,1)
                x = 10*log10(C.(m{1}).(sch{s,1})(groups{g,2}));
                plot(ax, sort(x), linspace(0,1,numel(x)), sch{s,4}, 'Color',sch{s,3}, 'LineWidth',1.5, 'DisplayName', sch{s,2});
            end
            xlabel(ax,'uplink SNR [dB]'); ylabel(ax,'CDF'); title(ax, groups{g,1});
            if g == 1, legend(ax,'Location','southeast','FontSize',8); end
        end
        title(tl, sprintf('SNR with 3 training entries (%s mode), noise %s', m{1}, strrep(nzName,'_','\_')));
        exportgraphics(f, fullfile(outDir, sprintf('snr_cdf_K3_%s.png', m{1})), 'Resolution',110);
        close(f);
    end
end

function write_readme(outDir, nz, T, modes, Ks, N, step, train, U, ch, specular, runtime)
    L = {sprintf('# Serving with AoA estimated from few EMS acquisitions (noise `%s`)', nz.name), '', ...
         sprintf('Generated %s by `experiments/Main_AoA_Serving.m` (branch aoa-serving), runtime %.0f min.', ...
                 char(datetime('now','Format','yyyy-MM-dd HH:mm')), runtime), '', ...
         '## What is simulated', ...
         sprintf('- **Serving codebook:** N = %d entries per EMS, phase step between adjacent elements from %g to %g deg/element in steps of %g (entry %d = zero phase = baseline). For lambda/2 spacing the beam is steered to sin(theta) = sin(theta_specular) + step/180; specular direction of EMS1 / EMS2 = %.1f / %.1f deg.', ...
                 N, step(1), step(end), step(2)-step(1), N/2, specular(1), specular(2)), ...
         sprintf('- **Training:** K entries spread evenly over the codebook (K = %s), e.g. K = 3 -> entries %s. `diag`: acquisition k sets both EMSs to training entry k (K acquisitions); `grid`: all pairs (K^2 acquisitions).', ...
                 mat2str(Ks), mat2str(train(3))), ...
         sprintf('- **AoA estimation:** one joint t-SNE chart (azimuth from EMS1 and EMS2) of all %d users, %d anchors (15%% of the users in front of both surfaces); one chart per training acquisition, then local fusion weighted by a label-free quality (no ground truth of the served users). Noise: %s', ...
                 U.n, ch.nA, nz.text), ...
         '- **Look-up tables (angle -> entry, per EMS):** `physics` = steering law with the known BS direction; `data` = majority best entry of the 7 anchors with the closest true angle (offline calibration with the labelled users).', ...
         '- **Link:** uplink, 23 dBm, noise -92 dBm, best rank-1 beamforming: SNR = P/N0 * lambda_max(H''H) with the TRUE channel H = direct + EMS1 + EMS2. SE = log2(1 + SNR). Coverage = share of users with SNR >= 0 dB (also 10 dB in the CSV).', ...
         '- **Compared schemes:** `baseline` (zero phase on both EMSs), `best_training` (best of the training configurations measured for the user, no angle estimation), `proposed_*` (t-SNE AoA + LUT), `genie_*` (true AoA + LUT: the cost of the angle error), `full_sweep` (best of all N^2 configurations, upper bound).', ...
         '- **Scored users:** unlabelled users only; groups: all, in front of both surfaces, without a direct path.', ''};
    for m = modes
        L{end+1} = sprintf('## Mode `%s`: median SNR gain vs baseline [dB] / mean SE [bit/s/Hz] / coverage SNR>=0 dB [%%]', m{1}); %#ok<AGROW>
        for g = unique(T.group,'stable').'
            L{end+1} = '';  L{end+1} = sprintf('**%s**', g{1}); %#ok<AGROW>
            L{end+1} = ['| scheme | ' strjoin(arrayfun(@(K) sprintf('K=%d', K), Ks, 'UniformOutput', false), ' | ') ' |']; %#ok<AGROW>
            L{end+1} = ['|---|' repmat('---|', 1, numel(Ks))]; %#ok<AGROW>
            for s = {'baseline','best_training','proposed_physics','proposed_data','genie_physics','genie_data','full_sweep'}
                c = cell(1,numel(Ks));
                for k = 1:numel(Ks)
                    q = T(strcmp(T.mode,m{1}) & strcmp(T.group,g{1}) & strcmp(T.scheme,s{1}) & T.K==Ks(k),:);
                    c{k} = sprintf('%.1f / %.2f / %.0f', q.dsnr_median_dB, q.se_mean, q.coverage_0dB);
                end
                L{end+1} = sprintf('| %s | %s |', s{1}, strjoin(c, ' | ')); %#ok<AGROW>
            end
        end
        q = T(strcmp(T.mode,m{1}) & strcmp(T.scheme,'baseline') & strcmp(T.group,'all users'),:);
        L{end+1} = ''; %#ok<AGROW>
        L{end+1} = sprintf('AoA error of the estimate (users in front of both, unlabelled), RMSE EMS1/EMS2 [deg]: %s', ...
            strjoin(arrayfun(@(i) sprintf('K=%d: %.1f/%.1f', q.K(i), q.aoa_rmse_ems1(i), q.aoa_rmse_ems2(i)), 1:height(q), 'UniformOutput', false), ', ')); %#ok<AGROW>
        L{end+1} = ''; %#ok<AGROW>
    end
    L = [L, {'Figures: `serving_vs_acquisitions_<mode>.png` (SNR gain, SE, coverage vs number of training acquisitions), `snr_cdf_K3_<mode>.png` (SNR CDFs with 3 training entries). All numbers: `summary.csv`.'}];
    fid = fopen(fullfile(outDir,'README.md'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
end

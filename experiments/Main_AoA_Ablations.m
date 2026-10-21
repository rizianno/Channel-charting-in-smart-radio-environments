%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AoA ABLATION STUDY
%
% Isolates, one factor at a time, what drives the AoA estimation (see the READMEs it writes):
%   users / charts : S1 joint (front of both) | S2 one 1-D chart per EMS | S3 joint, all users |
%                    S4 per-EMS on the S1 users (3x3, F0 only)
%   input feature  : F0 covariance per sweep + fusion | F1 stacked covariances |
%                    F2 sweep changes | F3 OFF subtraction
%   direct path    : on | off
%   codebook size  : 1x1 ... 5x5 per scenario (trend vs number of acquisitions)
%   codebook range : narrow (0..40 deg/element) | wide (-180..180 deg/element)
%   noise          : noiseless | estimated channels with 0 / 30 / 60 dB pilot gain
% Each t-SNE chart is computed once (F0 charts are shared by all codebook sizes) and cached in
% cache_ablations/ (git-ignored); outputs go to results_aoa/ablations/.
%
% Environment variables:
%   ABL_STAGE     calibration (default): S1, F0, direct on, 3x3, all four noise models
%                 matrix      : every scenario and codebook size, noise model ABL_NOISE
%                 evaluate    : redo only the post-processing of ABL_EVAL (calibration|<noise>)
%   ABL_NOISE     noiseless | est_G0dB | est_G30dB | est_G60dB (matrix stage)
%   ABL_WORKERS   parallel workers (default 16)
%   ABL_KS        restrict codebook sizes, e.g. "3" or "1,3"      (testing)
%   ABL_SCENARIOS restrict scenarios by name substrings, e.g. "S1-joint,S2-perEMS"
%   ABL_CODEBOOK  narrow | wide (default: both; charts already cached are not recomputed)
%   ABL_MAX_USERS use only the first N dataset users               (testing)
%   ABL_OUT       output root relative to the repo (default results_aoa/ablations)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();
params.target = 'aoa';
CAT    = abl_catalog();
stage  = envs('ABL_STAGE', 'calibration');
nW     = str2double(envs('ABL_WORKERS', '16'));
outRoot  = fullfile(root, envs('ABL_OUT', fullfile('results_aoa','ablations')));
cacheDir = fullfile(root, 'cache_ablations');
if ~exist(cacheDir,'dir'), mkdir(cacheDir); end

%% ------------------------------------------------------------------------
% DATA
%% ------------------------------------------------------------------------
Data = load_dataset(params);
maxU = str2double(envs('ABL_MAX_USERS', ''));
if ~isnan(maxU)
    for fld = {'R_xy','H_LOS_1bounce','H_I_RIS1','H_I_RIS2'}
        x = Data.(fld{1});  Data.(fld{1}) = x(1:maxU,:,:);
    end
end
U = abl_users(Data, params);
fprintf('Users: %d (front of both %d, EMS1 %d, EMS2 %d)\n', U.n, nnz(U.both), ...
        nnz(U.front(:,1)), nnz(U.front(:,2)));

%% ------------------------------------------------------------------------
% WHAT TO RUN
%% ------------------------------------------------------------------------
scns = CAT.scenarios;
switch stage
    case 'calibration'
        noises = CAT.noise;
        scns   = scns(strcmp({scns.structure},'S1') & strcmp({scns.feature},'F0') & [scns.direct] & strcmp({scns.codebook},'narrow'));
        scns.Ks = 3;
        stageDir = 'calibration';
    case 'matrix'
        noises = CAT.noise(strcmp({CAT.noise.name}, envs('ABL_NOISE', '')));
        assert(~isempty(noises), 'Set ABL_NOISE to one of: %s', strjoin({CAT.noise.name}, ', '));
        stageDir = noises.name;
    case 'evaluate'
        stageDir = envs('ABL_EVAL', 'calibration');
        if strcmp(stageDir,'calibration')
            noises = CAT.noise;
            scns   = scns(strcmp({scns.structure},'S1') & strcmp({scns.feature},'F0') & [scns.direct] & strcmp({scns.codebook},'narrow'));
            scns.Ks = 3;
        else
            noises = CAT.noise(strcmp({CAT.noise.name}, stageDir));
        end
    otherwise
        error('Unknown ABL_STAGE %s', stage);
end
ks = envs('ABL_KS', '');
if ~isempty(ks)
    kk = str2double(strsplit(ks, ','));
    for i = 1:numel(scns), scns(i).Ks = intersect(scns(i).Ks, kk); end
    scns = scns(arrayfun(@(s) ~isempty(s.Ks), scns));
end
cbf = envs('ABL_CODEBOOK', '');
if ~isempty(cbf), scns = scns(strcmp({scns.codebook}, cbf)); end
flt = envs('ABL_SCENARIOS', '');
if ~isempty(flt)
    pats = strtrim(strsplit(flt, ','));
    scns = scns(cellfun(@(n) any(cellfun(@(p) contains(n,p), pats)), {scns.name}));
end
fprintf('Stage %s: %d noise model(s), %d scenario(s)\n', stage, numel(noises), numel(scns));

%% ------------------------------------------------------------------------
% RUN (or load) THE t-SNE CHARTS
%% ------------------------------------------------------------------------
cacheFile = fullfile(cacheDir, sprintf('%s_U%d.mat', stageDir, U.n));
R = containers.Map();
if isfile(cacheFile)
    C0 = load(cacheFile);
    for i = 1:numel(C0.keys), R(C0.keys{i}) = C0.vals{i}; end
    fprintf('Loaded %d cached charts from %s\n', R.Count, cacheFile);
end

if ~strcmp(stage, 'evaluate')
    pool = gcp('nocreate');
    if isempty(pool), pool = parpool('Processes', nW); end
    for nz = noises
        % unique jobs of this noise model
        J = [];
        for s = scns
            charts = abl_charts(U, s.structure);
            for K = s.Ks, J = [J, abl_jobs(s, K, charts, CAT, nz.name)]; end %#ok<AGROW>
        end
        [~, iu] = unique({J.key}, 'stable');
        J = J(iu);
        J = J(~cellfun(@(k) isKey(R,k), {J.key}));
        fprintf('\n[%s] %d charts to compute (%.1f min elapsed)\n', nz.name, numel(J), toc(t0)/60);

        Fmap = containers.Map();
        getF = @(d, phi) get_features(Fmap, Data, params, d, phi, nz);
        futs = parallel.FevalFuture.empty;  fkeys = {};
        next = 1;  done = 0;
        while next <= numel(J) || ~isempty(futs)
            while next <= numel(J) && numel(futs) < 2*nW
                j = J(next);  next = next + 1;
                D = abl_job_input(j, getF);
                futs(end+1) = parfeval(pool, @abl_chart, 1, D, j.chart.labels, j.chart.nA, ...
                                       j.chart.valid, params.tsne_perplexity, params.K_fraction); %#ok<SAGROW>
                fkeys{end+1} = j.key; %#ok<SAGROW>
            end
            try
                [idx, out] = fetchNext(futs);
                R(fkeys{idx}) = out;
            catch ME
                idx = find(arrayfun(@(f) strcmp(f.State,'finished') && ~isempty(f.Error), futs), 1);
                fprintf(2, 'Job %s failed: %s\n', fkeys{idx}, ME.message);
            end
            futs(idx) = [];  fkeys(idx) = [];  done = done + 1;
            if mod(done, 25) == 0 || (next > numel(J) && isempty(futs))
                save_cache(cacheFile, R);
                fprintf('  %d/%d charts done (%.1f min)\n', done, numel(J), toc(t0)/60);
            end
        end
        clear Fmap getF
    end
end

%% ------------------------------------------------------------------------
% EVALUATION, README, TREND, COMPARISON
%% ------------------------------------------------------------------------
Tall = {};
for nz = noises
    for s = scns
        charts = abl_charts(U, s.structure);
        if strcmp(stageDir,'calibration')
            sDir = fullfile(outRoot, 'calibration', nz.name, s.name);
        else
            sDir = fullfile(outRoot, stageDir, s.name);
        end
        Ts = {};
        for K = s.Ks
            Ts{end+1} = abl_evaluate(s, K, U, charts, R, CAT, nz.name, fullfile(sDir, sprintf('K%d',K))); %#ok<SAGROW>
        end
        Ts = vertcat(Ts{:});
        abl_trend(Ts, s, sDir);
        abl_readme(s, CAT, nz, U, charts, scns, sDir, Ts);
        Ts.noise = repmat({nz.name}, height(Ts), 1);
        Tall{end+1} = Ts; %#ok<SAGROW>
        fprintf('evaluated %s / %s\n', nz.name, s.name);
    end
end
Tall = vertcat(Tall{:});
stageOut = fullfile(outRoot, stageDir);
if strcmp(stageDir,'calibration')
    abl_compare(Tall, stageOut, 'noise');
else
    abl_compare(Tall, stageOut, 'structure');
end
write_index(stageOut, stageDir, scns, noises, CAT, U, Tall, toc(t0)/60);
fprintf('\n=== DONE (%.1f min) -> %s ===\n', toc(t0)/60, stageOut);


%% ========================================================================
function v = envs(name, default)
    v = getenv(name);
    if isempty(v), v = default; end
end

function F = get_features(Fmap, Data, params, direct, phi, nz)
    k = abl_stkey(direct, phi);
    if ~isKey(Fmap, k)
        st = struct('direct', direct, 'phi', phi, 'key', k);
        Fmap(k) = abl_state_features(Data, params, st, nz);
    end
    F = Fmap(k);
end

function save_cache(file, R)
    keys = R.keys;  vals = R.values; %#ok<NASGU>
    save(file, 'keys', 'vals', '-v7.3');
end

function write_index(outDir, stageDir, scns, noises, CAT, U, T, runtime)
    if ~exist(outDir,'dir'), mkdir(outDir); end
    L = {};
    if strcmp(stageDir,'calibration')
        L{end+1} = '# AoA ablations: noise calibration';
        L{end+1} = '';
        L{end+1} = ['Reference scenario only (S1 joint chart, F0 covariance per sweep + fusion, ' ...
                    'direct path on, 3x3 codebook) under each noise model, to choose the noise ' ...
                    'model of the full ablation matrix.'];
    else
        L{end+1} = sprintf('# AoA ablations, noise model `%s`', stageDir);
        L{end+1} = '';
        L{end+1} = sprintf('Noise: %s', noises.text);
    end
    L{end+1} = '';
    L{end+1} = sprintf('Users: %d (in front of both surfaces: %d; EMS1: %d; EMS2: %d). Runtime %.0f min.', ...
                       U.n, nnz(U.both), nnz(U.front(:,1)), nnz(U.front(:,2)), runtime);
    L{end+1} = '';
    L{end+1} = '## Factors';
    for s = CAT.structures, L{end+1} = sprintf('- **%s** (`%s`): %s', s.id, s.tok, s.text); end %#ok<AGROW>
    for f = CAT.features,   L{end+1} = sprintf('- **%s** (`%s`): %s', f.id, f.tok, f.text); end %#ok<AGROW>
    L{end+1} = '- **Direct path** `dir-on` / `dir-off` (ideal removal).';
    for c = CAT.codebooks, L{end+1} = sprintf('- **Codebook range `%s`** (folder suffix `%s`): %s.', c.id, c.tok, c.text); end %#ok<AGROW>
    for n = noises, L{end+1} = sprintf('- **Noise `%s`**: %s', n.name, n.text); end %#ok<AGROW>
    L{end+1} = '';
    L{end+1} = '## Results at a glance (pooled EMS1+EMS2 median |error| [deg], subset `common`)';
    L{end+1} = '';
    L{end+1} = '| noise | scenario | method | K | acquisitions | pooled median | pooled mean | EMS1 median | EMS2 median |';
    L{end+1} = '|---|---|---|---|---|---|---|---|---|';
    M = T(ismember(T.method, {'local_orig','local_free','chart'}) & strcmp(T.subset,'common'), :);
    keys_ = unique(strcat(M.noise,'|',M.scenario,'|',M.method,'|',string(M.K)), 'stable');
    for k = keys_.'
        r = M(strcmp(strcat(M.noise,'|',M.scenario,'|',M.method,'|',string(M.K)), k), :);
        g = @(e) r(strcmp(r.ems,e), :);
        p = g('pooled'); e1 = g('EMS1'); e2 = g('EMS2');
        L{end+1} = sprintf('| %s | [%s](%s/) | %s | %d | %d | %.1f | %.1f | %.1f | %.1f |', p.noise{1}, ...
            p.scenario{1}, rel(stageDir, p.noise{1}, p.scenario{1}), p.method{1}, p.K, p.acquisitions, ...
            p.median_deg, p.mean_deg, e1.median_deg, e2.median_deg); %#ok<AGROW>
    end
    L{end+1} = '';
    L{end+1} = 'Each scenario folder has its own README (what changes and which scenarios isolate each factor), `trend.png` (error vs codebook size) and one `K<k>/` folder per codebook size. `comparison_common.png` / `comparison_natural.png` compare all scenarios.';
    fid = fopen(fullfile(outDir,'README.md'), 'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
end

function p = rel(stageDir, noise, scenario)
    if strcmp(stageDir,'calibration'), p = [noise '/' scenario]; else, p = scenario; end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PER-USER AoA ESTIMATES OF THE SEPARATE PER-EMS CHARTS, EVERY K (cached charts only, no channels)
%
% Each EMS trained alone (the other EMS absent), full-NLOS users it reaches; training informed or
% nested (coarse-to-fine order), K = 1..Kmax. Per user, the chart with the highest label-free
% trustworthiness. Same caches, training lists, anchors and selection as
% Main_AoA_Serving_Deployments (its per_user.csv has K = 1, 3, 6 only; this file has every K).
% Scored users: reached by the EMS and not anchors.
% Environment: EAP_SIZES ("60r,60,120,180,240"), EAP_KMAX (15), EAP_FOOT (0.15), EAP_NOISE (est_G60dB)
% Outputs: results_aoa/serving_deployments/<noise>/per_user_aoa.csv
%          .../training_lists.csv: training codeword of slot k (and its direction u = sin theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
t0 = tic;

params = init_parameters();  params.target = 'aoa';
CAT   = abl_catalog();
nz    = CAT.noise(strcmp({CAT.noise.name}, envs('EAP_NOISE','est_G60dB')));
sizes = strtrim(strsplit(envs('EAP_SIZES','60r,60,120,180,240'), ','));
Kmax  = str2double(envs('EAP_KMAX','15'));
foot  = str2double(envs('EAP_FOOT','0.15'));

Data = load_dataset(params);
U    = abl_users(Data, params);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);  specular = -azBS;
popU = find(~U.has_direct);  nU = numel(popU);                  % full NLOS, as the deployment run
reach = [any(Data.H_I_RIS1(popU,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(popU,:,:) ~= 0, [2 3])] & U.front(popU,:);
clear Data
ch = struct('users',{}, 'nA',{}, 'lab',{}, 'uA',{});
isAnch = false(nU,2);
for e = 1:2
    m = find(reach(:,e));  nA = round(params.Supervision_percentage * numel(m));
    ch(e).users = m;  ch(e).nA = nA;  ch(e).lab = U.az(popU(m), e);  ch(e).uA = sind(ch(e).lab(1:nA));
    isAnch(m(1:nA), e) = true;
end
azTrue = U.az(popU, :);

rows = {};  tl = {};
for sz = sizes
    real_ = endsWith(sz{1}, 'r');  M = str2double(erase(sz{1}, 'r'));
    tag = sprintf('M%d%s', M, repmat('_raytraced', 1, real_));
    wu = 0.886 * 2 / M;  N = ceil(2/wu);  N = N + mod(N,2);
    step = ((1:N) - N/2) * 360/N;
    for e = 1:2
        ub = mod(sind(specular(e)) + step/180 + 1, 2) - 1;
        no = nest_order(N);
        lists = {'nested', no(1:Kmax); ...
                 'informed', informed_order(ch(e).uA, ub, foot, Kmax, max(2, round(0.02*ch(e).nA)))};
        cf = fullfile(root,'cache_ablations', sprintf('single_%s_e%d_fullnlos_%s.mat', nz.name, e, tag));
        Ra = load_cache(cf);
        for t = 1:size(lists,1)
            li = lists{t,2};
            tl = [tl; repmat({tag, sprintf('EMS%d', e), lists{t,1}}, Kmax, 1), ...
                  num2cell([(1:Kmax).', li(:), ub(li).', repmat([N, sind(specular(e))], Kmax, 1)])]; %#ok<AGROW>
            for K = 1:Kmax
                if ~all(arrayfun(@(i) isKey(Ra, sprintf('%d', i)), li(1:K)))
                    fprintf('[%s] EMS%d %s: K = %d not cached, stop\n', tag, e, lists{t,1}, K);  break
                end
                outs = Ra(sprintf('%d', li(1)));
                for q = 2:K, outs(q) = Ra(sprintf('%d', li(q))); end
                est = select_chart([outs.est], [outs.tw], ch(e).nA, ch(e).lab);
                ix = ch(e).users;
                err = abs(mod(est - azTrue(ix,e) + 180, 360) - 180);
                rows = [rows; repmat({tag, sprintf('EMS%d', e), lists{t,1}, K}, numel(ix), 1), ...
                        num2cell([popU(ix), U.xy(popU(ix),:), isAnch(ix,e), ~isAnch(ix,e), azTrue(ix,e), est, err])]; %#ok<AGROW>
            end
        end
        fprintf('[%s] EMS%d done (%.1f min)\n', tag, e, toc(t0)/60);
    end
end
outDir = fullfile(root,'results_aoa','serving_deployments', nz.name);
writetable(cell2table(rows, 'VariableNames', {'surface','view','training','K','user','x_m','y_m','anchor','scored', ...
    'az_true_deg','az_est_deg','aoa_err_deg'}), fullfile(outDir,'per_user_aoa.csv'));
writetable(cell2table(tl, 'VariableNames', {'surface','view','training','k','codeword','u','N','sin_spec'}), ...
           fullfile(outDir,'training_lists.csv'));
fprintf('\n=== DONE (%.1f min) ===\n', toc(t0)/60);

%% ========================================================================
function v = envs(name, default)
    v = getenv(name);  if isempty(v), v = default; end
end

function x = select_chart(E, QT, nA, lab)
% per user, the chart with the highest label-free trustworthiness (anchors keep their labels)
    [n, K] = size(E);  [~, k] = max(QT, [], 2);  x = E(sub2ind([n K], (1:n).', k));  x(1:nA) = lab(1:nA);
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

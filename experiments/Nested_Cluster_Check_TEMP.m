%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEMPORARY DIAGNOSTIC: why does EMS2's AoA estimate (nested training, per-user selection)
% suddenly improve at K = 5 on the real 60 x 60 surface?
%   1) beam direction of each nested training entry, per EMS (steering law)
%   2) EMS2 median AoA error per spatial user cluster, for K = 1..10 (selection), from the
%      cached nested charts (no new t-SNE), and which chart each cluster's users select
% Output: printed + results_aoa/serving_nlos_size_nestedFusion_TEMP/cluster_check.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
params = init_parameters();  params.target = 'aoa';
Data = load_dataset(params);  U = abl_users(Data, params);
idx = find(~U.has_direct);  nU = numel(idx);
reachE = [any(Data.H_I_RIS1(idx,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(idx,:,:) ~= 0, [2 3])] & U.front(idx,:);
azBS = compute_ems_azimuth(U.BS_position(1:2), U.EMS_position, params);  spec = -azBS;
M = 60;  N = ceil(2/(0.886*2/M));  N = N + mod(N,2);  step = ((1:N) - N/2) * 360/N;
f = [1/2 1];  for lev = 2:12, f = [f, (1:2:2^lev-1)/2^lev]; if numel(unique(round(f*N))) >= N, break; end, end %#ok<AGROW>
ord = unique(min(max(round(f*N),1),N), 'stable');  ord = ord(1:10);
out = {};
out{end+1} = sprintf('Nested entries (real 60x60, N = %d): %s', N, mat2str(ord));
for e = 1:2
    s = mod(sind(spec(e)) + step(ord)/180 + 1, 2) - 1;          % sin of the beam direction (wrapped)
    out{end+1} = sprintf('EMS%d beam azimuth per entry: %s deg', e, mat2str(round(asind(s)))); %#ok<SAGROW>
end
% EMS2 chart users and their spatial clusters (as in the cluster summary)
e = 2;  m = find(reachE(:,e));  nA = round(params.Supervision_percentage * numel(m));
xy = U.xy(idx(m),:);  az = U.az(idx(m),e);  lab = az;
cl = dbscan(xy, 4, 8);
C = load(fullfile(root,'cache_ablations','nested_est_G60dB_M60_raytraced.mat'));  R = containers.Map(C.keys, C.vals);
unl = (nA+1:numel(m)).';
big = find(accumarray(cl(cl>0), 1) >= 40);
for k = big.'
    s = unl(cl(unl) == k);
    out{end+1} = sprintf('cluster %d: %d users (scored), centre (%.0f, %.0f) m, azimuth %.0f..%.0f deg', ...
                         k, numel(s), mean(xy(s,1)), mean(xy(s,2)), min(az(s)), max(az(s))); %#ok<SAGROW>
end
for K = 1:10
    E = [];  Q = [];
    for i = ord(1:K)
        o = R(sprintf('est_G60dB|M60_raytraced|nest|e%d|%d', e, i));  E = [E o.est]; Q = [Q o.tw]; %#ok<AGROW>
    end
    [~, kk] = max(Q, [], 2);  sel = E(sub2ind(size(E), (1:size(E,1)).', kk));
    err = abs(mod(sel - lab + 180, 360) - 180);
    line = sprintf('K=%2d (adds entry %2d):', K, ord(K));
    for k = big.'
        s = unl(cl(unl) == k);
        cnt = histcounts(kk(s), 0.5:1:K+0.5);  [~, dom] = max(cnt);
        line = [line sprintf('  cl%d %5.1f deg (picks entry %d for %3.0f%%)', k, median(err(s)), ord(dom), 100*max(cnt)/numel(s))]; %#ok<AGROW>
    end
    out{end+1} = line; %#ok<SAGROW>
end
fprintf('%s\n', out{:});
fid = fopen(fullfile(root,'results_aoa','serving_nlos_size_nestedFusion_TEMP','cluster_check.txt'),'w');
fprintf(fid, '%s\n', out{:});  fclose(fid);

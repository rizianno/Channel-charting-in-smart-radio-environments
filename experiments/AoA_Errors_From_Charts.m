%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AoA ERRORS OF THE PER-EMS t-SNE CHARTS OF Main_AoA_Serving_NLoS_Size (from the cached charts)
% For every population (full NLOS, banana = no LoS), surface and number of training sweeps K:
% the fused azimuth estimates (label-free local fusion, as in the serving experiment) and their
% median / RMSE of |error| for EMS1, EMS2 and joint sqrt(e1^2 + e2^2).
%   EMS e : users of chart e that are not its anchors
%   joint : users in both charts and anchor in neither
% No new t-SNE: reads cache_ablations/nlos_size_<noise>_<tag>.mat.
% Output: results_aoa/serving_nlos_size_sameEntryTraining/aoa_errors.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
params = init_parameters();  params.target = 'aoa';
nzName = 'est_G60dB';
Data = load_dataset(params);
U = abl_users(Data, params);
Nb = params.Par.Tx.ArrSize;  Nr = size(Data.H_LOS_1bounce,3);
Pd = squeeze(sum(sum(abs(Data.H_LOS_1bounce).^2,2),3)) / (Nb*Nr);
dd = sqrt(sum(([U.xy, params.UE_height*ones(U.n,1)] - U.BS_position).^2, 2));
los = U.has_direct & 10*log10(Pd ./ ((3e8/30e9)./(4*pi*dd)).^2) > -36;

rows = {};
for pop = {'fullnlos','nolos'}
    if strcmp(pop{1},'nolos'), idx = find(~los); else, idx = find(~U.has_direct); end
    nU = numel(idx);  isNone = ~U.has_direct(idx);
    reachE = [any(Data.H_I_RIS1(idx,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(idx,:,:) ~= 0, [2 3])] & U.front(idx,:);
    ch = struct('users',{}, 'nA',{}, 'labels',{});  isAnch = false(nU,2);
    for e = 1:2
        m = find(reachE(:,e));  nA = round(params.Supervision_percentage * numel(m));
        ch(e).users = m;  ch(e).nA = nA;  ch(e).labels = U.az(idx(m), e);
        isAnch(m(1:nA), e) = true;
    end
    for tagc = {'M60_raytraced','M60','M120','M180','M240'}
        tag = tagc{1};  M = sscanf(tag, 'M%d');
        ctag = tag;  if strcmp(pop{1},'nolos'), ctag = ['nolos_' tag]; end
        f = fullfile(root,'cache_ablations', sprintf('nlos_size_%s_%s.mat', nzName, ctag));
        if ~isfile(f), continue; end
        C = load(f);  R = containers.Map(C.keys, C.vals);
        N = ceil(2/(0.886*2/M));  N = N + mod(N,2);
        spread = @(s) unique(round((2*(1:s) - 1) * N / (2*s)));
        for K = 1:8
            tr = spread(K);  az = nan(nU,2);
            for e = 1:2
                outs = R(sprintf('%s|%s|e%d|%d', nzName, ctag, e, tr(1)));
                for s = 2:numel(tr), outs(s) = R(sprintf('%s|%s|e%d|%d', nzName, ctag, e, tr(s))); end
                E3 = cat(3, outs.est);
                if numel(tr) == 1, est = E3; else, [~, est] = abl_fuse(E3, [outs.q_free], ch(e).nA, ch(e).labels); end
                az(ch(e).users, e) = est;
            end
            err = mod(az - U.az(idx,:) + 180, 360) - 180;
            for sub = {'full NLOS', isNone; 'all users of the population', true(nU,1)}.'
                if strcmp(pop{1},'fullnlos') && strcmp(sub{1},'all users of the population'), continue; end
                s1 = reachE(:,1) & ~isAnch(:,1) & sub{2};  s2 = reachE(:,2) & ~isAnch(:,2) & sub{2};
                sj = all(reachE,2) & ~any(isAnch,2) & sub{2};
                j  = sqrt(err(sj,1).^2 + err(sj,2).^2);
                st = @(x) [median(abs(x)), sqrt(mean(x.^2))];
                rows(end+1,:) = {pop{1}, sub{1}, tag, K, nnz(s1), nnz(s2), nnz(sj), ...
                                 st(err(s1,1)), st(err(s2,2)), st(j)}; %#ok<SAGROW>
                % same, split: users reachable by both EMSs vs by this EMS only
                o1 = s1 & ~reachE(:,2);  o2 = s2 & ~reachE(:,1);
                rows(end+1,:) = {pop{1}, [sub{1} ' | reachable by both'], tag, K, nnz(sj), nnz(sj), nnz(sj), ...
                                 st(err(sj,1)), st(err(sj,2)), st(j)}; %#ok<SAGROW>
                rows(end+1,:) = {pop{1}, [sub{1} ' | reachable by this EMS only'], tag, K, nnz(o1), nnz(o2), 0, ...
                                 st(err(o1,1)), st(err(o2,2)), [NaN NaN]}; %#ok<SAGROW>
            end
        end
    end
end
rows = [rows(:,1:7), num2cell(cell2mat(rows(:,8:10)))];
T = cell2table(rows, 'VariableNames', {'population','users','surface','K','n_ems1','n_ems2','n_joint', ...
    'ems1_median_deg','ems1_rmse_deg','ems2_median_deg','ems2_rmse_deg','joint_median_deg','joint_rmse_deg'});
writetable(T, fullfile(root,'results_aoa','serving_nlos_size_sameEntryTraining','aoa_errors.csv'));
disp(T);

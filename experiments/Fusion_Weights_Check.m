%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DIAGNOSTIC: why does fusing several t-SNE charts give a worse AoA than one chart?
% For the full-NLOS per-EMS charts of Main_AoA_Serving_NLoS_Size (cached, no new t-SNE), for
% each training set K: every single chart's median AoA error, its fusion weight (global and the
% spread of the local weights), and the error of the fused estimate. Same-entry training sets.
% Output: results_aoa/serving_nlos_size_sameEntryTraining/fusion_weights_check.csv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
params = init_parameters();  params.target = 'aoa';
nzName = 'est_G60dB';
Data = load_dataset(params);
U = abl_users(Data, params);
idx = find(~U.has_direct);  nU = numel(idx);
reachE = [any(Data.H_I_RIS1(idx,:,:) ~= 0, [2 3]), any(Data.H_I_RIS2(idx,:,:) ~= 0, [2 3])] & U.front(idx,:);
rows = {};
for tagc = {'M60_raytraced','M240'}
    tag = tagc{1};  M = sscanf(tag, 'M%d');
    C = load(fullfile(root,'cache_ablations', sprintf('nlos_size_%s_%s.mat', nzName, tag)));
    R = containers.Map(C.keys, C.vals);
    N = ceil(2/(0.886*2/M));  N = N + mod(N,2);
    spread = @(s) unique(round((2*(1:s) - 1) * N / (2*s)));
    for e = 1:2
        m = find(reachE(:,e));  nA = round(params.Supervision_percentage * numel(m));
        lab = U.az(idx(m), e);  unl = nA+1:numel(m);
        for K = 1:8
            tr = spread(K);
            outs = R(sprintf('%s|%s|e%d|%d', nzName, tag, e, tr(1)));
            for s = 2:numel(tr), outs(s) = R(sprintf('%s|%s|e%d|%d', nzName, tag, e, tr(s))); end
            E3 = squeeze(cat(3, outs.est));  if isvector(E3), E3 = E3(:); end      % n x K
            err = abs(mod(E3 - lab + 180, 360) - 180);
            medc = median(err(unl,:), 1);                                       % per-chart median error
            Q = [outs.q_free];
            gq = arrayfun(@(c) prctile(Q(unl,c), 90), 1:numel(tr));  wg = gq / sum(gq);
            Wl = max(Q,1e-6);  Wl = Wl ./ sum(Wl,2);
            if numel(tr) > 1
                [g, l] = abl_fuse(reshape(E3, [], 1, numel(tr)), Q, nA, lab);
            else
                g = E3;  l = E3;
            end
            eg = median(abs(mod(g(unl) - lab(unl) + 180, 360) - 180));
            el = median(abs(mod(l(unl) - lab(unl) + 180, 360) - 180));
            rows(end+1,:) = {tag, e, K, mat2str(tr), mat2str(round(medc,1)), mat2str(round(wg,3)), ...
                             sprintf('%.3f..%.3f', min(Wl(:)), max(Wl(:))), min(medc), mean(medc), eg, el}; %#ok<SAGROW>
        end
    end
end
T = cell2table(rows, 'VariableNames', {'surface','ems','K','training_entries','single_chart_median_err_deg', ...
    'global_weights','local_weight_range','best_single_deg','typical_single_deg','global_fused_deg','local_fused_deg'});
writetable(T, fullfile(root,'results_aoa','serving_nlos_size_sameEntryTraining','fusion_weights_check.csv'));
disp(T(:, {'surface','ems','K','single_chart_median_err_deg','global_weights','best_single_deg','typical_single_deg','global_fused_deg','local_fused_deg'}));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST: AoA FROM ESTIMATED POSITIONS vs DIRECT AoA ESTIMATION
%
% Takes the position charts of Main_Position_CodebookSize (cached in
% results_position_codebook_charts.mat, no new t-SNE), fuses them per codebook size as
% that script does, converts the estimated positions into the azimuth of the user seen
% from EMS1 and EMS2, and scores these angles. They are compared with charts that use
% the azimuths directly as t-SNE labels (branch aoa-estimation / aoa-ablations):
%   AOA_DIRECT_REF : results.mat of the direct-AoA run with the SAME settings as the
%                    position run (noiseless, 1e-6 regularization, narrow 3x3), i.e.
%                    results_aoa/coarse_3x3/results.mat on branch aoa-estimation
%   AOA_ABL_DIR    : results_aoa/ablations/est_G60dB of branch aoa-ablations
%                    (direct AoA at K = 1..5, estimated channels with 60 dB pilot gain:
%                    different noise model, secondary comparison)
% Scored users: in front of both surfaces and anchor in none of the compared runs.
% Azimuth: from each surface normal, positive toward +y (codebook-gradient axis), as
% compute_ems_azimuth on the AoA branches; errors wrapped to [-180, 180).
% Output: results_position_codebook/position_to_aoa/ (README.md, summary.csv)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(root,'src')));
params = init_parameters();
Data   = load_dataset(params);
state  = prepare_data(Data, params);
Len = state.Len;  nA = params.supervised_Len;

EMSpos = [Data.RIS1_position; Data.RIS2_position];
normal = [-1 0 0; 1 0 0];          % outward normals (EMS1 faces -x, EMS2 faces +x)
axis_m = [ 0 1 0; 0 1 0];          % codebook-gradient axis
azim  = @(xy) [atan2d((xy - EMSpos(1,1:2))*axis_m(1,1:2).', (xy - EMSpos(1,1:2))*normal(1,1:2).'), ...
               atan2d((xy - EMSpos(2,1:2))*axis_m(2,1:2).', (xy - EMSpos(2,1:2))*normal(2,1:2).')];
wrap  = @(x) mod(x + 180, 360) - 180;
azT   = azim(state.R_xy);
front = [(state.R_xy - EMSpos(1,1:2))*normal(1,1:2).' > 0, (state.R_xy - EMSpos(2,1:2))*normal(2,1:2).' > 0];
both  = all(front, 2);

% anchors: position run = first nA users; AoA runs = first 15% of the users in front of both
idxB = find(both);  aoaAnch = false(Len,1);  aoaAnch(idxB(1:round(0.15*numel(idxB)))) = true;
posAnch = false(Len,1);  posAnch(1:nA) = true;
S = both & ~posAnch & ~aoaAnch;
fprintf('Scored users: %d\n', nnz(S));

%% position charts -> fusion -> azimuth
C0 = load(fullfile(root,'results_position_codebook_charts.mat'));
wide = @(K) sort(mod((0:K-1)*360/K + 180, 360) - 180);
LEV.narrow = {0, [0 40], [0 20 40], [0 40/3 80/3 40], [0 10 20 30 40]};
LEV.wide   = arrayfun(wide, 1:5, 'UniformOutput', false);
unl = nA+1:Len;
rows = {};
for r = {'narrow','wide'}
    for K = 1:5
        [a,b] = ndgrid(LEV.(r{1}){K}, LEV.(r{1}){K});
        [~, loc] = ismember(round([a(:) b(:)],6), C0.cfgs, 'rows');
        R = [C0.CH{loc}];
        rp.PosEstimates        = cat(3, R.xy);
        rp.ErrorPerCodebook    = arrayfun(@(o) mean(o.err(unl)), R).';
        rp.CT_PerUser          = [R.tc1];
        rp.TW_PerUser          = [R.tc2];
        rp.CT_PerCodebook_mean = mean(rp.CT_PerUser(unl,:)).';
        rp.CT_PerCodebook_p90  = prctile(rp.CT_PerUser(unl,:),90).';
        rp.TW_PerCodebook_mean = mean(rp.TW_PerUser(unl,:)).';
        rp.TW_PerCodebook_p90  = prctile(rp.TW_PerUser(unl,:),90).';
        rl = local_weighting(rp, state, params);
        rg = global_weighting(rp, state, params);
        rows(end+1,:) = [{sprintf('position -> AoA, %s', r{1}), K, 'local_fusion'}, num2cell(score(wrap(azim(rl.xy) - azT), S))]; %#ok<SAGROW>
        rows(end+1,:) = [{sprintf('position -> AoA, %s', r{1}), K, 'global_fusion'}, num2cell(score(wrap(azim(rg.xy) - azT), S))]; %#ok<SAGROW>
        sc = zeros(numel(loc), 6);
        for c = 1:numel(loc), sc(c,:) = score(wrap(azim(rp.PosEstimates(:,:,c)) - azT), S); end
        rows(end+1,:) = [{sprintf('position -> AoA, %s', r{1}), K, 'typical_single'}, num2cell(mean(sc,1))]; %#ok<SAGROW>
    end
end

%% direct AoA, same settings (noiseless, 1e-6, narrow 3x3)
ref = getenv('AOA_DIRECT_REF');
if ~isempty(ref) && isfile(ref)
    A = load(ref, 'est', 'Y');
    assert(isequal(size(A.Y,1), numel(idxB)), 'direct AoA run: unexpected user set');
    for m = {'local','global'}
        E = nan(Len,2);  E(idxB,:) = A.est.(m{1});
        assert(max(abs(wrap(A.Y - azT(idxB,:))), [], 'all') < 1e-6, 'azimuth definitions differ');
        rows(end+1,:) = [{'direct AoA, narrow (old run, same settings)', 3, [m{1} '_fusion']}, ...
                         num2cell(score(wrap(E - azT), S))]; %#ok<SAGROW>
    end
end

%% direct AoA, ablation runs (estimated channels, 60 dB pilot gain)
abl = getenv('AOA_ABL_DIR');
scn = {'S1-joint_F0-cov-per-sweep_dir-on', 'S3-joint-allUE_F0-cov-per-sweep_dir-on'};
lbl = {'direct AoA (est. channels 60 dB), users in front of both', 'direct AoA (est. channels 60 dB), all users in t-SNE'};
if ~isempty(abl) && isfolder(abl)
    for s = 1:numel(scn)
        for cb = {'', '_cb-wide'}
            for K = 1:5
                f = fullfile(abl, [scn{s} cb{1}], sprintf('K%d',K), 'per_user_local_orig.csv');
                if ~isfile(f), continue; end
                P = readtable(f);
                E = [P.az1_est_deg P.az2_est_deg];
                rng_ = 'narrow'; if ~isempty(cb{1}), rng_ = 'wide'; end
                rows(end+1,:) = [{sprintf('%s, %s', lbl{s}, rng_), K, 'local_fusion'}, ...
                                 num2cell(score(wrap(E - azT), S))]; %#ok<SAGROW>
            end
        end
    end
end

T = cell2table(rows, 'VariableNames', {'method','K','fusion','EMS1_rmse_deg','EMS2_rmse_deg', ...
               'joint_rmse_deg','EMS1_median_deg','EMS2_median_deg','joint_median_deg'});
outDir = fullfile(root,'results_position_codebook','position_to_aoa');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir,'summary.csv'));
disp(T);

%% README
L = {'# AoA from estimated positions vs direct AoA estimation (test)', '', ...
     sprintf(['The position charts of `Main_Position_CodebookSize` (no new t-SNE) are fused per codebook size, ' ...
              'the estimated positions are converted into the azimuth seen from EMS1 and EMS2, and scored on ' ...
              '%d users (in front of both surfaces, anchor in none of the compared runs). Generated %s by ' ...
              '`experiments/Position_to_AoA.m`.'], nnz(S), char(datetime('now','Format','yyyy-MM-dd HH:mm'))), '', ...
     '**Fair comparison** (same noiseless setting, 1e-6 regularization, same 3x3 narrow codebook): `position -> AoA, narrow` K=3 vs `direct AoA, narrow (old run, same settings)`. The ablation rows use estimated channels (60 dB pilot gain), so they are a secondary comparison.', '', ...
     '| method | K | fusion | EMS1 RMSE | EMS2 RMSE | joint RMSE | EMS1 median | EMS2 median | joint median |', ...
     '|---|---|---|---|---|---|---|---|---|'};
for i = 1:height(T)
    L{end+1} = sprintf('| %s | %dx%d | %s | %.1f | %.1f | %.1f | %.1f | %.1f | %.1f |', T.method{i}, T.K(i), T.K(i), ...
        T.fusion{i}, T{i,4:9}); %#ok<SAGROW>
end
L{end+1} = '';
L{end+1} = 'Angles in degrees. Joint = sqrt(e1^2 + e2^2) per user. Position run: all 3200 users, 480 anchors (first users in dataset order); direct AoA runs: users in front of both surfaces (or all users), 305 anchors.';
fid = fopen(fullfile(outDir,'README.md'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
fprintf('\n=== DONE -> %s ===\n', outDir);

function s = score(E, S)
    e1 = E(S,1);  e2 = E(S,2);  j = sqrt(e1.^2 + e2.^2);
    s = [sqrt(mean(e1.^2)), sqrt(mean(e2.^2)), sqrt(mean(j.^2)), median(abs(e1)), median(abs(e2)), median(j)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TEST: AoA FROM ESTIMATED POSITIONS vs DIRECT AoA ESTIMATION
%
% Takes the position charts of Main_Position_CodebookSize (cached in
% results_position_codebook_charts.mat, no new t-SNE), fuses them per codebook size as
% that script does, converts the estimated positions into the azimuth of the user seen
% from EMS1 and EMS2 ("position proxy"), and compares these angles with charts that use
% the azimuths directly as t-SNE labels (branch aoa-ablations), K = 1..5, both ranges:
%   AOA_ABL_ROOT/noiseless_1e6 : direct AoA with the SAME setting as the position run
%                                (noiseless, old 1e-6 regularization)  -> fair comparison
%   AOA_ABL_ROOT/est_G60dB     : direct AoA with estimated channels (60 dB pilot gain)
% Direct AoA runs used: per-sweep charts + fusion with the direct path (as the position
% run), (a) users in front of both surfaces, (b) all users in t-SNE.
% Scored users: in front of both surfaces and anchor in none of the compared runs.
% Azimuth: from each surface normal, positive toward +y (codebook-gradient axis), as
% compute_ems_azimuth on the AoA branches; errors wrapped to [-180, 180).
% Output: results_position_codebook/position_to_aoa/ (README.md, summary.csv, figures)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
set(0,'DefaultFigureVisible','off');
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

% anchors: position run = first nA users; AoA joint charts = first 15% of the users in front of both
idxB = find(both);  aoaAnch = false(Len,1);  aoaAnch(idxB(1:round(0.15*numel(idxB)))) = true;
posAnch = false(Len,1);  posAnch(1:nA) = true;
S = both & ~posAnch & ~aoaAnch;
fprintf('Scored users: %d\n', nnz(S));

rows = {};
add = @(rows, src, users, rng_, K, meth, E) [rows; ...
      {src, users, rng_, K, K^2, meth, 'joint'},  num2cell(stat(sqrt(E(S,1).^2 + E(S,2).^2))); ...
      {src, users, rng_, K, K^2, meth, 'pooled'}, num2cell(stat(abs([E(S,1); E(S,2)])))];

%% position proxy: position charts -> fusion -> azimuth
C0 = load(fullfile(root,'results_position_codebook_charts.mat'));
wide = @(K) sort(mod((0:K-1)*360/K + 180, 360) - 180);
LEV.narrow = {0, [0 40], [0 20 40], [0 40/3 80/3 40], [0 10 20 30 40]};
LEV.wide   = arrayfun(wide, 1:5, 'UniformOutput', false);
unl = nA+1:Len;
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
        src = 'position proxy (noiseless, 1e-6)';
        rows = add(rows, src, 'all users', r{1}, K, 'local_fusion',  wrap(azim(rl.xy) - azT));
        rows = add(rows, src, 'all users', r{1}, K, 'global_fusion', wrap(azim(rg.xy) - azT));
        for c = 1:numel(loc)
            rows = add(rows, src, 'all users', r{1}, K, sprintf('single_%d',c), ...
                       wrap(azim(rp.PosEstimates(:,:,c)) - azT));
        end
    end
end

%% direct AoA (ablation runs): local fusion of per-sweep charts with the direct path
abl = getenv('AOA_ABL_ROOT');
scn = {'S1-joint_F0-cov-per-sweep_dir-on', 'users in front of both'; ...
       'S3-joint-allUE_F0-cov-per-sweep_dir-on', 'all users'};
nzs = {'noiseless_1e6', 'direct AoA (noiseless, 1e-6)'; 'est_G60dB', 'direct AoA (est. channels, 60 dB)'};
for n = 1:size(nzs,1)
    for s = 1:size(scn,1)
        for cb = {'', '_cb-wide'}
            for K = 1:5
                f = fullfile(abl, nzs{n,1}, [scn{s,1} cb{1}], sprintf('K%d',K), 'per_user_local_orig.csv');
                if ~isfile(f), continue; end
                P = readtable(f);
                assert(max(abs(wrap([P.az1_true_deg P.az2_true_deg] - azT)), [], 'all') < 1e-6, ...
                       'azimuth definitions differ');
                rng_ = 'narrow'; if ~isempty(cb{1}), rng_ = 'wide'; end
                rows = add(rows, nzs{n,2}, scn{s,2}, rng_, K, 'local_fusion', ...
                           wrap([P.az1_est_deg P.az2_est_deg] - azT));
            end
        end
    end
end

T = cell2table(rows, 'VariableNames', {'source','users','range','K','acquisitions','method','error', ...
               'median_deg','mean_deg','p90_deg','rmse_deg'});
% typical single configuration of the position proxy = mean over the K^2 single charts
sg = startsWith(T.method,'single_');
G  = groupsummary(T(sg,:), {'source','users','range','K','acquisitions','error'}, 'mean', ...
                  {'median_deg','mean_deg','p90_deg','rmse_deg'});
G  = renamevars(G, strcat('mean_', {'median_deg','mean_deg','p90_deg','rmse_deg'}), ...
                {'median_deg','mean_deg','p90_deg','rmse_deg'});
G.method = repmat({'typical_single'}, height(G), 1);
T = [T(~sg,:); G(:, T.Properties.VariableNames)];

outDir = fullfile(root,'results_position_codebook','position_to_aoa');
if ~exist(outDir,'dir'), mkdir(outDir); end
writetable(T, fullfile(outDir,'summary.csv'));

%% figures
L1 = {'position proxy (noiseless, 1e-6)','all users',[0 0 0]; ...
      'direct AoA (noiseless, 1e-6)','users in front of both',[0 0.45 0.74]; ...
      'direct AoA (noiseless, 1e-6)','all users',[0.30 0.75 0.93]};
L2 = {'direct AoA (noiseless, 1e-6)','users in front of both',[0 0.45 0.74]; ...
      'direct AoA (noiseless, 1e-6)','all users',[0.30 0.75 0.93]; ...
      'direct AoA (est. channels, 60 dB)','users in front of both',[0.85 0.33 0.10]; ...
      'direct AoA (est. channels, 60 dB)','all users',[0.93 0.69 0.13]};
fig(T, L1, 'AoA from estimated positions vs direct AoA estimation (same noiseless setting, local fusion)', ...
    fullfile(outDir,'proxy_vs_direct.png'), S);
fig(T, L2, 'Direct AoA estimation: noiseless vs estimated channels (60 dB pilot gain), local fusion', ...
    fullfile(outDir,'direct_noiseless_vs_noisy.png'), S);

%% README
L = {'# AoA from estimated positions vs direct AoA estimation (test)', '', ...
     sprintf(['The position charts of `Main_Position_CodebookSize` (no new t-SNE) are fused per codebook size, ' ...
              'the estimated positions are converted into the azimuth seen from EMS1 and EMS2 ("position proxy"), ' ...
              'and compared with direct AoA charts (azimuths as t-SNE labels, per-sweep charts + local fusion, ' ...
              'direct path kept), on %d users (in front of both surfaces, anchor in none of the runs). ' ...
              'Generated %s by `experiments/Position_to_AoA.m`.'], nnz(S), char(datetime('now','Format','yyyy-MM-dd HH:mm'))), '', ...
     '- `proxy_vs_direct.png`: **fair comparison**, all noiseless with the 1e-6 regularization (the position run''s setting).', ...
     '- `direct_noiseless_vs_noisy.png`: direct AoA, noiseless vs estimated channels (60 dB pilot gain).', ...
     '- Rows of each figure: joint error sqrt(e1^2+e2^2) and per-surface error (EMS1 and EMS2 pooled); columns: median, mean, 90th percentile, RMSE; solid = narrow phase steps (0..40 deg/element), dashed = wide (-180..180).', '', ...
     '## Joint error, local fusion [deg]: median / mean / p90 / RMSE', '', ...
     '| source | users | range | 1x1 | 2x2 | 3x3 | 4x4 | 5x5 |', '|---|---|---|---|---|---|---|---|'};
for l = unique(strcat(T.source,'|',T.users,'|',T.range), 'stable').'
    p = strsplit(l{1}, '|');
    cells = cell(1,5);
    for K = 1:5
        q = T(strcmp(T.source,p{1}) & strcmp(T.users,p{2}) & strcmp(T.range,p{3}) & T.K==K & ...
              strcmp(T.method,'local_fusion') & strcmp(T.error,'joint'), :);
        if isempty(q), cells{K} = '-'; else, cells{K} = sprintf('%.1f / %.1f / %.1f / %.1f', q{1,8:11}); end
    end
    L{end+1} = sprintf('| %s | %s | %s | %s |', p{1}, p{2}, p{3}, strjoin(cells, ' | ')); %#ok<SAGROW>
end
L = [L, {'', 'All statistics (also global fusion and the typical single configuration of the position proxy, joint and pooled errors) are in `summary.csv`.'}];
fid = fopen(fullfile(outDir,'README.md'),'w'); fprintf(fid, '%s\n', strjoin(L, newline)); fclose(fid);
fprintf('\n=== DONE -> %s ===\n', outDir);


function s = stat(x)
    s = [median(x), mean(x), prctile(x,90), sqrt(mean(x.^2))];
end

function fig(T, lines_, ttl, file, S)
    M = {'median_deg','median'; 'mean_deg','mean'; 'p90_deg','90th percentile'; 'rmse_deg','RMSE'};
    E = {'joint','joint error \surd(e_1^2+e_2^2)'; 'pooled','per-surface error (EMS1, EMS2 pooled)'};
    f = figure('Position',[50 50 1800 820]);
    tl = tiledlayout(f, 2, 4, 'TileSpacing','compact');
    for r = 1:2
        for m = 1:4
            ax = nexttile(tl); hold(ax,'on'); grid(ax,'on'); box(ax,'on');
            for l = 1:size(lines_,1)
                for rg = {'narrow','wide'}
                    q = sortrows(T(strcmp(T.source,lines_{l,1}) & strcmp(T.users,lines_{l,2}) & ...
                         strcmp(T.range,rg{1}) & strcmp(T.method,'local_fusion') & strcmp(T.error,E{r,1}),:), 'K');
                    if isempty(q), continue; end
                    ls = '-o'; if strcmp(rg{1},'wide'), ls = '--s'; end
                    plot(ax, q.acquisitions, q.(M{m,1}), ls, 'Color',lines_{l,3}, 'LineWidth',1.5, 'MarkerSize',4, ...
                         'DisplayName', sprintf('%s, %s, %s', lines_{l,1}, lines_{l,2}, rg{1}));
                end
            end
            set(ax,'XScale','log'); xticks(ax,[1 4 9 16 25]); xlim(ax,[0.9 28]);
            xlabel(ax,'EMS configurations (K^2)'); ylabel(ax, sprintf('%s [deg]', M{m,2}));
            title(ax, sprintf('%s: %s', M{m,2}, E{r,1}));
            if r == 1 && m == 4, legend(ax,'Location','eastoutside','FontSize',7); end
        end
    end
    title(tl, sprintf('%s, %d users', ttl, nnz(S)));
    exportgraphics(f, file, 'Resolution', 110); close(f);
end

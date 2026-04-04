%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Main_DoubleEMS_Trajectory.m
%
%  Generates:
%    Fig. A  — Trajectory preview on the ground-truth grid + speed/accel
%              kinematics plots
%    Fig. B  — Side-by-side trajectory in channel-chart space:
%              EMS-assisted vs Direct-only, with per-point error coloring
%
%  The user picks a codeword in one of two ways:
%    CODEBOOK_SOURCE = 'manual'  →  codebook_idx is used directly
%    CODEBOOK_SOURCE = 'auto'    →  best-index is loaded from the saved
%                                   data of Main_Double_EMS_BestCodebook.m
%                                   (field = AUTO_METRIC + '_' + AUTO_AGG)
%
%  Run-or-load logic
%  -----------------
%  If a previously saved embedding file is found the expensive computation
%  (covariance + dissimilarity + t-SNE) is skipped entirely and the script
%  jumps straight to trajectory generation and figure plotting.
%  Only the codebook index and trajectory seed must agree with the saved
%  file (a mismatch triggers a warning and recomputes).
%
%  Dissimilarity : same as chosen in the BestCodebook run (or override here)
%  Embedding     : semi-supervised t-SNE  (evaluate_channel_configuration)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% ========================================================================
%%% SECTION 0 — STANDARD INIT  (mirrors Main_Multi_Slot_Fusion)
%%% ========================================================================

if isempty(mfilename('fullpath'))
    root = pwd;
else
    root = fileparts(mfilename('fullpath'));
end

addpath(genpath(fullfile(root, '..', 'src')));
params = init_parameters();

root   = get_project_root();
Data   = load_dataset(params);
state  = prepare_data(Data, params);

%% ========================================================================
%%% SECTION 1 — EXPERIMENT-SPECIFIC CHOICES
%%% ========================================================================

% -------------------------------------------------------------------------
%  Codebook selection
%  CODEBOOK_SOURCE = 'manual' : use codebook_idx below directly
%  CODEBOOK_SOURCE = 'auto'   : load best index from BestCodebook saved data
% -------------------------------------------------------------------------
% CODEBOOK_SOURCE = 'auto';    % 'manual' | 'auto'
CODEBOOK_SOURCE = 'manual';    % 'manual' | 'auto'
codebook_idx    = 61;        % used only when CODEBOOK_SOURCE = 'manual'

%  When CODEBOOK_SOURCE = 'auto': which field of bestIdxStruct to use?
%    AUTO_METRIC : 'CT' | 'TW' | 'ER'
%    AUTO_AGG    : 'mean' | 'median' | 'pX'
AUTO_METRIC = 'CT';
AUTO_AGG    = 'pX';

% -------------------------------------------------------------------------
%  BestCodebook reference parameters
%  (only needed when CODEBOOK_SOURCE = 'auto'; must match the BestCodebook
%   run that produced the saved .mat file we want to load from)
% -------------------------------------------------------------------------
bc_dissim_type      = 'LE';    % dissim type used in BestCodebook run
bc_supPct           = 0.15;    % supervision percentage used there
bc_nCodes           = [];      % leave [] to auto-detect from codebook

% -------------------------------------------------------------------------
%  Dissimilarity metric for THIS run's embeddings
% -------------------------------------------------------------------------
dissim_type = 'LE';            % 'LE' | 'CMD'

% -------------------------------------------------------------------------
%  Supervision
% -------------------------------------------------------------------------
Supervision_percentage = 0.15;
supervised_Len = round(params.Total_samples_used * Supervision_percentage);

% -------------------------------------------------------------------------
%  Trajectory parameters  (matching original main_physical_vs_latent.m)
% -------------------------------------------------------------------------
totalTime  = 5000;   % [s]    total simulation time
dt         = 5;      % [s]    time step
v_max      = [];     % [m/s]  leave [] → auto-selected from grid spacing
a_max      = 0.2;    % [m/s²] small → smoother path
startIdx   = 1;      % starting node index in R_xy
K_fallback = 8;      % K-NN fallback edges for graph connectivity

% Random seed for the trajectory walk.
%   'default' resets to MATLAB factory state (Mersenne Twister, seed 0).
%   The original code had rng('shuffle') commented out and no prior randn
%   calls before the trajectory, so the effective seed was 'default'.
%   Set to a positive integer for a reproducible non-default path, or []
%   to leave the generator in whatever state it currently has.
rng_seed   = 'default';   % 'default' | positive integer | []

% -------------------------------------------------------------------------
%  Run-or-load override
%  Set FORCE_RECOMPUTE = true to rerun channel charting even if a saved
%  embedding file already exists (useful after changing parameters).
% -------------------------------------------------------------------------
FORCE_RECOMPUTE = false;

% -------------------------------------------------------------------------
%  Visualisation parameters
% -------------------------------------------------------------------------
ERR_THR = 10;   % [m]  error threshold for coloring trajectory points

% -------------------------------------------------------------------------
%  Derived variables from state / params
% -------------------------------------------------------------------------
Len     = state.Len;
R_xy    = state.R_xy;
ArrSize = params.Par.Tx.ArrSize;
N_real  = params.N_real;
M       = params.M;
N       = params.N;

% Codebook (for index lookup and phi recovery)
codebook     = generate_codebook(params);
nCodes       = codebook.nCodebooks;
codebook_all = codebook.all;    % [nCodes x 2]

% Dissimilarity function handle
if strcmpi(dissim_type, 'LE')
    dissimFn = @CompDismilarityLogEuclidean;
else
    dissimFn = @CompDismilarity_CMD;
end

%% ========================================================================
%% SECTION 2 — OUTPUT FOLDERS
%% ========================================================================

samplesStr = sprintf('%d_Samples', params.Total_samples_used);

figFolder = fullfile(root, 'Figures', samplesStr, ...
            ['DoubleEMS_Trajectory_' dissim_type], filesep);
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

%% ========================================================================
%% SECTION 3 — RESOLVE CODEBOOK INDEX
%% ========================================================================

if strcmpi(CODEBOOK_SOURCE, 'manual')
    % -----------------------------------------------------------------
    %  Manual selection: use the index supplied in Section 1
    % -----------------------------------------------------------------
    codebook_idx_used = codebook_idx;
    fprintf('=== Codebook source: MANUAL ===\n');
    fprintf('Using codebook index %d\n\n', codebook_idx_used);

else
    % -----------------------------------------------------------------
    %  Auto selection: load bestIdxStruct from BestCodebook saved data
    % -----------------------------------------------------------------
    fprintf('=== Codebook source: AUTO (%s_%s) ===\n', AUTO_METRIC, AUTO_AGG);

    if isempty(bc_nCodes)
        bc_nCodes = nCodes;   % assume same codebook size
    end

    bc_fileTag  = sprintf('DoubleEMS_%s_samp%d_super%.2f_nCodes%d', ...
                          bc_dissim_type, params.Total_samples_used, ...
                          bc_supPct, bc_nCodes);
    bc_dataFile = fullfile(root, 'merged_results_doubleEMS', bc_fileTag, ...
                           ['MergedData_' bc_fileTag '.mat']);

    assert(isfile(bc_dataFile), ...
        ['BestCodebook data file not found:\n  %s\n' ...
         'Run Main_Double_EMS_BestCodebook.m first.'], bc_dataFile);

    loaded      = load(bc_dataFile, 'bestIdxStruct', 'codebook_all');
    fieldName   = [AUTO_METRIC '_' AUTO_AGG];
    assert(isfield(loaded.bestIdxStruct, fieldName), ...
        'Field ''%s'' not found in bestIdxStruct. Check AUTO_METRIC / AUTO_AGG.', fieldName);

    codebook_idx_used = loaded.bestIdxStruct.(fieldName);
    codebook_all      = loaded.codebook_all;   % use the saved codebook

    fprintf('Loaded best index %d  (field: %s)\n\n', codebook_idx_used, fieldName);
end

phi1_used = codebook_all(codebook_idx_used, 1);
phi2_used = codebook_all(codebook_idx_used, 2);
fprintf('phi1 = %.4f rad  (%.1f deg)   phi2 = %.4f rad  (%.1f deg)\n\n', ...
        phi1_used, rad2deg(phi1_used), phi2_used, rad2deg(phi2_used));

%% ========================================================================
%%% SECTION 4 — CHANNEL CHARTING: RUN OR LOAD
%%% ========================================================================

fileTag  = sprintf('Traj_DoubleEMS_%s_samp%d_super%.2f_cidx%d', ...
                   dissim_type, params.Total_samples_used, ...
                   Supervision_percentage, codebook_idx_used);
dataFile = fullfile(figFolder, [fileTag '.mat']);

if isfile(dataFile) && ~FORCE_RECOMPUTE
    % ------------------------------------------------------------------
    %  Fast path: load saved embeddings
    % ------------------------------------------------------------------
    fprintf('=== Saved embeddings found ===\nLoading from:\n  %s\n\n', dataFile);
    tmp = load(dataFile, 'xy_tsne_EMS', 'xy_tsne_DIR', 'codebook_idx_saved');

    % Warn if saved index doesn't match
    if tmp.codebook_idx_saved ~= codebook_idx_used
        warning(['Saved embeddings were computed for codebook index %d, ' ...
                 'but codebook_idx_used = %d. Recomputing.'], ...
                 tmp.codebook_idx_saved, codebook_idx_used);
        clear tmp;
        embeddings_loaded = false;
    else
        xy_tsne_EMS       = tmp.xy_tsne_EMS;
        xy_tsne_DIR       = tmp.xy_tsne_DIR;
        embeddings_loaded = true;
        fprintf('Embeddings loaded. Skipping t-SNE computation.\n\n');
    end
else
    if FORCE_RECOMPUTE
        fprintf('FORCE_RECOMPUTE=true: ignoring any cached embeddings.\n\n');
    end
    embeddings_loaded = false;
end

if ~embeddings_loaded
    % ------------------------------------------------------------------
    %  Full computation: build H_tot, covariance, dissimilarity, t-SNE
    % ------------------------------------------------------------------
    fprintf('=== Computing channel charting embeddings ===\n\n');

    % --- EMS-assisted channel -----------------------------------------
    fprintf('[1/2] EMS-assisted channel (phi1=%.3f, phi2=%.3f) ...\n', ...
            phi1_used, phi2_used);

    phase1 = kron(phi1_used*(1:M), ones(1,N));
    phase2 = kron(phi2_used*(1:M), ones(1,N));
    Theta1 = diag(exp(-1j * phase1));
    Theta2 = diag(exp(-1j * phase2));

    H_ems1 = zeros(Len, ArrSize, size(state.H_i_EMS1, 3));
    H_ems2 = zeros(Len, ArrSize, size(state.H_i_EMS2, 3));
    for uu = 1:Len
        H_ems1(uu,:,:) = state.H_o_EMS1 * Theta1 * squeeze(state.H_i_EMS1(uu,:,:));
        H_ems2(uu,:,:) = state.H_o_EMS2 * Theta2 * squeeze(state.H_i_EMS2(uu,:,:));
    end
    H_tot = state.H_dir + H_ems1 + H_ems2;

    [~, ~, ~, ~, ~, xy_tsne_EMS] = evaluate_channel_configuration( ...
        H_tot, R_xy, supervised_Len, ArrSize, N_real, dissimFn);

    % --- Direct-only (No EMS) baseline --------------------------------
    fprintf('[2/2] Direct-only baseline ...\n');
    [~, ~, ~, ~, ~, xy_tsne_DIR] = evaluate_channel_configuration( ...
        state.H_dir, R_xy, supervised_Len, ArrSize, N_real, dissimFn);

    % --- Save embeddings ----------------------------------------------
    codebook_idx_saved = codebook_idx_used;
    save(dataFile, ...
         'xy_tsne_EMS', 'xy_tsne_DIR', ...
         'codebook_idx_saved', 'phi1_used', 'phi2_used', ...
         'supervised_Len', 'Supervision_percentage', ...
         'dissim_type', 'params', '-v7.3');
    fprintf('\nEmbeddings saved to:\n  %s\n\n', dataFile);
end

%% ========================================================================
%% SECTION 5 — TRAJECTORY GENERATION
%% ========================================================================

% Auto-select v_max from grid spacing if not specified
if isempty(v_max)
    nearDist = min(pdist2(R_xy, R_xy) + diag(inf(size(R_xy,1), 1)), [], 2);
    v_max    = 1.2 * median(nearDist) / dt;
    fprintf('Auto-selected  v_max = %.2f m/s  (d_step = %.2f m)\n', ...
            v_max, v_max*dt);
end

% Set random seed BEFORE printing so a string 'default' doesn't crash fprintf
if ~isempty(rng_seed)
    rng(rng_seed);    % 'default' restores MATLAB factory state (seed 0)
end

if ischar(rng_seed)
    fprintf('=== Generating trajectory (seed=''%s'', T=%.0f s, dt=%.1f s) ===\n', ...
            rng_seed, totalTime, dt);
else
    fprintf('=== Generating trajectory (seed=%d, T=%.0f s, dt=%.1f s) ===\n', ...
            rng_seed, totalTime, dt);
end

[trajIdx, trajXY, speed, accel, timeVec] = ...
    buildTrajectory_smooth_local(R_xy, totalTime, dt, v_max, a_max, startIdx, K_fallback);

fprintf('Trajectory built: %d steps (%.1f s)\n\n', numel(trajIdx), timeVec(end));

%% ========================================================================
%%% SECTION 6 — POSITION ESTIMATION ALONG TRAJECTORY
%%%   tsne_sem aligns the latent chart to physical coordinates via the
%%%   supervised anchors → xy_tsne values are directly in metres and can
%%%   be compared to the physical trajectory.
%%% ========================================================================

estXY_EMS = xy_tsne_EMS(trajIdx, :);   % estimated physical positions (EMS)
estXY_DIR = xy_tsne_DIR(trajIdx, :);   % estimated physical positions (No EMS)

% Per-step positioning error (Euclidean distance in physical space)
err_EMS = sqrt(sum((estXY_EMS - trajXY).^2, 2));
err_DIR = sqrt(sum((estXY_DIR - trajXY).^2, 2));

% Good / bad masks
good_EMS = err_EMS <= ERR_THR;   bad_EMS = ~good_EMS;
good_DIR = err_DIR <= ERR_THR;   bad_DIR = ~good_DIR;

%% ========================================================================
%%% SECTION 7 — FIGURES
%%% ========================================================================
% % % % % % % fprintf('=== Generating figures ===\n');
% % % % % % % 
% % % % % % % % =========================================================================
% % % % % % % %  FIGURE A — Physical trajectory preview (map + directional arrows)
% % % % % % % % =========================================================================
% % % % % % % figA = figure('Position', [70 70 700 580]);
% % % % % % % hold on; box on; grid on;
% % % % % % % 
% % % % % % % scatter(R_xy(:,1), R_xy(:,2), 6, [0.83 0.83 0.83], 'filled');
% % % % % % % plot(trajXY(:,1), trajXY(:,2), '-', 'LineWidth', 1.8, ...
% % % % % % %      'Color', [0 0.45 0.74]);
% % % % % % % 
% % % % % % % % Arrows every ~4 s
% % % % % % % stepA = max(round(4 / dt), 1);
% % % % % % % idxA  = 1:stepA:(numel(timeVec)-1);
% % % % % % % u = trajXY(idxA+1, 1) - trajXY(idxA, 1);
% % % % % % % v = trajXY(idxA+1, 2) - trajXY(idxA, 2);
% % % % % % % quiver(trajXY(idxA,1), trajXY(idxA,2), u, v, 0, ...
% % % % % % %        'Color',       [0 0.45 0.74], ...
% % % % % % %        'LineWidth',   1, ...
% % % % % % %        'MaxHeadSize', 0.35, ...
% % % % % % %        'AutoScale',   'off');
% % % % % % % 
% % % % % % % xlabel('$X$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
% % % % % % % ylabel('$Y$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
% % % % % % % title(sprintf('Candidate trajectory  ($T$=%.0f s,  $v_{\\max}$=%.2f,  $a_{\\max}$=%.2f)', ...
% % % % % % %               totalTime, v_max, a_max), ...
% % % % % % %       'Interpreter', 'latex', 'FontSize', 14);
% % % % % % % axis equal; set(gca, 'FontSize', 14);
% % % % % % % 
% % % % % % % sgtitle('\textbf{Preview of physical trajectory (check before localisation)}', ...
% % % % % % %         'Interpreter', 'latex', 'FontSize', 14);
% % % % % % % 
% % % % % % % savefig(fullfile(figFolder, ['FigA_TrajPreview_' fileTag '.fig']));
% % % % % % % saveas(figA,    fullfile(figFolder, ['FigA_TrajPreview_' fileTag '.png']));

% =========================================================================
%  FIGURE B — Physical-space comparison: EMS-assisted vs Direct-only
%  Estimated positions plotted on the physical grid alongside the actual
%  trajectory.  Colour: orange/green = good (≤ ERR_THR m), red = bad.
% =========================================================================
figB = figure('Position', [60 60 1320 520]);

% --- Left panel: EMS-assisted --------------------------------------------
subplot(1, 2, 1); hold on; box on; grid on;

scatter(R_xy(:,1), R_xy(:,2), 6, [0.85 0.85 0.85], 'filled');
plot(trajXY(:,1), trajXY(:,2), '-', 'LineWidth', 2, 'Color', [0 0.45 0.74]);

scatter(estXY_EMS(good_EMS,1), estXY_EMS(good_EMS,2), ...
        18, 'filled', ...
        'MarkerFaceColor', [0.85 0.33 0.10], 'MarkerEdgeColor', 'k');
scatter(estXY_EMS(bad_EMS,1),  estXY_EMS(bad_EMS,2), ...
        28, 'filled', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'k');

axis equal;
xlabel('$X$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$Y$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
title({'EMS + direct channel'; ...
       sprintf('($\\varphi_1 = %.4f$ rad,  $\\varphi_2 = %.4f$ rad)', ...
               phi1_used, phi2_used)}, ...
      'Interpreter', 'latex', 'FontSize', 14);
legend({'Grid', 'Actual path', ...
        sprintf('Estimate ($\\leq %d$ m)', ERR_THR), ...
        sprintf('Estimate ($> %d$ m)',     ERR_THR)}, ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
set(gca, 'FontSize', 14);

% --- Right panel: Direct only --------------------------------------------
subplot(1, 2, 2); hold on; box on; grid on;

scatter(R_xy(:,1), R_xy(:,2), 6, [0.85 0.85 0.85], 'filled');
plot(trajXY(:,1), trajXY(:,2), '-', 'LineWidth', 2, 'Color', [0 0.45 0.74]);

scatter(estXY_DIR(good_DIR,1), estXY_DIR(good_DIR,2), ...
        18, 'filled', ...
        'MarkerFaceColor', [0.3 0.7 0.3], 'MarkerEdgeColor', 'k');
scatter(estXY_DIR(bad_DIR,1),  estXY_DIR(bad_DIR,2), ...
        28, 'filled', ...
        'MarkerFaceColor', [1 0 0], 'MarkerEdgeColor', 'k');

axis equal;
xlabel('$X$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
ylabel('$Y$ [m]', 'Interpreter', 'latex', 'FontSize', 16);
title('Direct channel only', 'Interpreter', 'latex', 'FontSize', 14);
legend({'Grid', 'Actual path', ...
        sprintf('Estimate ($\\leq %d$ m)', ERR_THR), ...
        sprintf('Estimate ($> %d$ m)',     ERR_THR)}, ...
       'Location', 'best', 'Interpreter', 'latex', 'FontSize', 13);
set(gca, 'FontSize', 14);

sgtitle(sprintf('\\bf Physical vs latent trajectories   (error threshold = %d m)', ERR_THR), ...
        'Interpreter', 'latex', 'FontSize', 14);

savefig(fullfile(figFolder, ['FigB_PhysVsLatent_' fileTag '.fig']));
saveas(figB,    fullfile(figFolder, ['FigB_PhysVsLatent_' fileTag '.png']));

fprintf('Figures saved to:\n  %s\n\n', figFolder);
fprintf('=== Done ===\n');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  LOCAL COPY of buildTrajectory_smooth
%%%  -----------------------------------------------------------------------
%%%  Kept here as a local function so this script is self-contained and is
%%%  guaranteed to use the exact original algorithm regardless of what is on
%%% the MATLAB path.  MATLAB always prefers a local function over a path
%%%  function when called from within the same file.
%%%
%%%  Algorithm is identical to src/buildTrajectory_smooth.m and to the
%%%  original main_physical_vs_latent.m local function.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [trajIdx, trajXY, speed, accel, timeVec] = ...
         buildTrajectory_smooth_local(R_xy, totalTime, dt, v_max, a_max, startIdx, K)
%     rng('shuffle')          % uncomment for a fresh random seed each run
    if nargin < 7, K = 4; end

    nTicks   = floor(totalTime / dt) + 1;
    d_max    = v_max * dt;
    a_max_dt = a_max * dt;          % max speed change per tick [m/s]

    %% radius graph + K-NN fallback
    D   = squareform(pdist(R_xy));
    adj = (D > 0) & (D <= d_max);  % radius graph
    if any(sum(adj, 2) == 0)        % connectivity fix
        [~, ord] = sort(D, 2);
        for i = 1:size(D, 1)
            adj(i, ord(i, 2:K+1)) = true;
        end
    end

    %% allocate
    trajIdx    = zeros(nTicks, 1);  trajIdx(1) = startIdx;
    speed      = zeros(nTicks, 1);
    accel      = zeros(nTicks, 1);
    prevDir    = [];                % previous heading unit-vector

    %% main loop
    for t = 2:nTicks
        cur  = trajIdx(t-1);
        nbrs = find(adj(cur, :));   % feasible neighbours by speed

        if isempty(nbrs)            % should not happen with K-NN fallback
            [~, nbrs] = min(D(cur, :));
        end

        v_prev    = (t > 2) * speed(t-1);   % 0 on first move
        bestScore = -inf;
        chosen    = [];

        % pick neighbour that respects accel & maximises heading alignment
        for n = nbrs(randperm(numel(nbrs)))
            dist = D(cur, n);   v = dist / dt;
            if abs(v - v_prev) > a_max_dt + 1e-12, continue; end
            if isempty(prevDir)
                score = 1;      % first segment → no heading preference
            else
                dirVec = (R_xy(n, :) - R_xy(cur, :)) / dist;
                score  = dirVec * prevDir';
            end
            if score > bestScore
                bestScore = score;
                chosen    = [n, dist, v];
            end
        end

        if isempty(chosen)          % relax acceleration constraint
            n    = nbrs(1);
            dist = D(cur, n);
            v    = dist / dt;
        else
            n    = chosen(1);
            dist = chosen(2);
            v    = chosen(3);
        end

        trajIdx(t) = n;
        speed(t)   = v;
        accel(t)   = (v - v_prev) / dt;
        prevDir    = (R_xy(n, :) - R_xy(cur, :)) / dist;
    end

    trajXY  = R_xy(trajIdx, :);
    timeVec = (0:nTicks-1).' * dt;
end

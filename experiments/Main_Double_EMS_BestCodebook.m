%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Main_DoubleEMS_BestCodebook.m
%
%  Generates:
%    Fig. 7  — CDF of −Continuity and −Trustworthiness (t-SNE, Double-EMS)
%              highlighting the specular codeword and the best codeword
%    Fig. 8  — CDF of Positioning Error (t-SNE, Double-EMS)
%              highlighting the specular codeword and the best codeword
%
%  For each metric the figure shows:
%    • Shaded envelope   — range across all static codewords
%    • Thin grey dashed  — every individual codeword
%    • Specular (φ=0,0)  — reference codeword
%    • Best codeword     — selected by 90th-percentile of the metric
%    • No EMS baseline   — direct channel only
%
%  Dissimilarity : 'LE' (Log-Euclidean) or 'CMD' (Correlation Matrix Dist.)
%  Embedding     : semi-supervised t-SNE
%
%  Run-or-load logic
%  -----------------
%  If a previously saved data file is found in outFolder the heavy
%  computation (codebook loop + baseline) is skipped entirely and the
%  script jumps straight to the CDF computation and figure generation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% ========================================================================
%% SECTION 0 — STANDARD INIT  (mirrors Main_Multi_Slot_Fusion)
%% ========================================================================

if isempty(mfilename('fullpath'))
    root = pwd;
else
    root = fileparts(mfilename('fullpath'));
end

addpath(genpath(fullfile(root, '..', 'src')));
params = init_parameters();

root          = get_project_root();
Data          = load_dataset(params);
state         = prepare_data(Data, params);

%% ========================================================================
%% SECTION 1 — EXPERIMENT-SPECIFIC CHOICES  (override shared params here)
%% ========================================================================

% --- Dissimilarity metric: 'LE' (Log-Euclidean) or 'CMD' ----------------
dissim_type = 'LE';

% --- Supervision --------------------------------------------------------
%     (overrides the 15% default in init_parameters)
Supervision_percentage = 0.15;
supervised_Len = round(params.Total_samples_used * Supervision_percentage);

% --- Derived convenience variables from state / params ------------------
Len      = state.Len;
R_xy     = state.R_xy;
ArrSize  = params.Par.Tx.ArrSize;
N_real   = params.N_real;
M        = params.M;
N        = params.N;

% --- Best-codeword selection -------------------------------------------
%     BEST_AGG controls which aggregation is used to highlight the "best"
%     codeword in the CDF figures AND to populate bestIdxStruct (used by
%     Main_DoubleEMS_Trajectory to auto-select a codebook).
%       'mean'   — mean of the per-user metric vector
%       'median' — median
%       'pX'     — bestPercentile-th percentile  (set below)
BEST_AGG       = 'pX';   % 'mean' | 'median' | 'pX'
bestPercentile = 90;     % only used when BEST_AGG = 'pX'

% --- Codebook  (phi_options = deg2rad(0:4:40), 11 x 11 = 121 entries) ---
codebook     = generate_codebook(params);
nCodes       = codebook.nCodebooks;
codebook_all = codebook.all;        % [nCodes x 2]
phi_options_1 = params.phi_options;
phi_options_2 = params.phi_options;

% --- Pick the dissimilarity function handle ----------------------------
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
            ['DoubleEMS_' dissim_type], filesep);
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

fileTag = sprintf('DoubleEMS_%s_samp%d_super%.2f_nCodes%d', ...
                  dissim_type, params.Total_samples_used, ...
                  Supervision_percentage, nCodes);

outFolder = fullfile(root, 'merged_results_doubleEMS', fileTag, filesep);
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

%% ========================================================================
%% SECTION 3 — RUN OR LOAD
%% ========================================================================

dataFile = fullfile(outFolder, ['MergedData_' fileTag '.mat']);

if isfile(dataFile)
    % ------------------------------------------------------------------
    %  Fast path: previously saved results found — skip all computation
    % ------------------------------------------------------------------
    fprintf('=== Saved data found ===\nLoading from:\n  %s\n\n', dataFile);
    load(dataFile);
    fprintf('Load complete. Skipping computation.\n\n');

else
    % ------------------------------------------------------------------
    %  Full computation path
    % ------------------------------------------------------------------
    fprintf('=== No saved data found. Running full computation. ===\n\n');

    %% ----------------------------------------------------------------
    %% 3A — STATIC CODEBOOK LOOP
    %% ----------------------------------------------------------------
    fprintf('=== Static codebook loop (%d entries) ===\n', nCodes);

    CTnegCell = cell(nCodes, 1);
    TWnegCell = cell(nCodes, 1);
    ERRcell   = cell(nCodes, 1);

    for idx = 1:nCodes
        phi1 = codebook_all(idx, 1);
        phi2 = codebook_all(idx, 2);
        fprintf('[%3d/%3d]  phi1 = %.3f rad  phi2 = %.3f rad\n', ...
                idx, nCodes, phi1, phi2);

        % Phase profiles and diagonal phase matrices
        phase1 = kron(phi1*(1:M), ones(1,N));
        phase2 = kron(phi2*(1:M), ones(1,N));
        Theta1 = diag(exp(-1j * phase1));
        Theta2 = diag(exp(-1j * phase2));

        % Build EMS-reflected channels using state fields
        H_ems1 = zeros(Len, ArrSize, size(state.H_i_EMS1, 3));
        H_ems2 = zeros(Len, ArrSize, size(state.H_i_EMS2, 3));
        for uu = 1:Len
            H_ems1(uu,:,:) = state.H_o_EMS1 * Theta1 * squeeze(state.H_i_EMS1(uu,:,:));
            H_ems2(uu,:,:) = state.H_o_EMS2 * Theta2 * squeeze(state.H_i_EMS2(uu,:,:));
        end

        % NOTE: no sqrt(8)/scaling factor — matches original script convention
        H_tot = state.H_dir + H_ems1 + H_ems2;

        % Evaluate CT, TW, error  (gain and dvec_dB discarded here)
        [CTn, TWn, ERR, ~, ~] = evaluate_channel_configuration( ...
            H_tot, R_xy, supervised_Len, ArrSize, N_real, dissimFn);

        CTnegCell{idx} = CTn;
        TWnegCell{idx} = TWn;
        ERRcell{idx}   = ERR;
    end
    fprintf('Static codebook loop complete.\n\n');

    %% ----------------------------------------------------------------
    %% 3B — NO-EMS BASELINE  (direct channel only, no scaling)
    %% ----------------------------------------------------------------
    fprintf('=== Computing No-EMS baseline ===\n');

    [CTdir, TWdir, err_dir, ~, ~] = evaluate_channel_configuration( ...
        state.H_dir, R_xy, supervised_Len, ArrSize, N_real, dissimFn);

    fprintf('No-EMS baseline done.\n\n');

    %% ----------------------------------------------------------------
    %% 3C — SAVE
    %% ----------------------------------------------------------------
    fprintf('Saving data to:\n  %s\n', dataFile);
    save(dataFile, ...
        'CTnegCell', 'TWnegCell', 'ERRcell', ...
        'CTdir', 'TWdir', 'err_dir', ...
        'codebook_all', 'phi_options_1', 'phi_options_2', ...
        'nCodes', 'supervised_Len', ...
        'params', 'dissim_type', 'Supervision_percentage', ...
        '-v7.3');
    fprintf('Data saved successfully.\n\n');
end

%% ========================================================================
%% SECTION 4 — CDF COMPUTATION
%%   (always runs — fast; kept outside the save/load block so bin edges
%%    can be adjusted freely without rerunning the simulation)
%% ========================================================================
fprintf('=== Computing empirical CDFs ===\n');

edgesCT = linspace(-1, 0, 201);
edgesER = linspace(0, 150, 201);

% Build [nCodes x 200] CDF matrices (one row per codeword)
% Note: MATLAB does not allow {: } indexing directly on a function call
% result, so cellfun and vertcat are split into two steps.
tmp_CT = cellfun(@(v) histcounts(v, edgesCT, 'Norm','cdf'), CTnegCell, 'uni', 0);
tmp_TW = cellfun(@(v) histcounts(v, edgesCT, 'Norm','cdf'), TWnegCell, 'uni', 0);
tmp_ER = cellfun(@(v) histcounts(v, edgesER, 'Norm','cdf'), ERRcell,   'uni', 0);
cdf_CT = vertcat(tmp_CT{:});
cdf_TW = vertcat(tmp_TW{:});
cdf_ER = vertcat(tmp_ER{:});

% Baseline CDFs (single row vectors)
cdf_CT_dir = histcounts(CTdir,    edgesCT, 'Norm','cdf');
cdf_TW_dir = histcounts(TWdir,    edgesCT, 'Norm','cdf');
cdf_ER_dir = histcounts(err_dir,  edgesER, 'Norm','cdf');

x_CT = edgesCT(2:end);
x_TW = x_CT;
x_ER = edgesER(2:end);

fprintf('CDFs computed.\n\n');

%% ========================================================================
%% SECTION 5 — BEST & SPECULAR CODEWORD INDICES
%%   All three aggregation methods are computed and stored.
%%   BEST_AGG (Section 1) selects which one is highlighted in the figures.
%%   CT/TW: negated, so HIGHER (closer to 0) = better  →  maximize
%%   ER:    raw error, so LOWER                = better  →  minimize
%% ========================================================================

% --- Compute aggregate vectors (one value per codeword) -----------------
CT_mean_vec   = cellfun(@mean,   CTnegCell);
TW_mean_vec   = cellfun(@mean,   TWnegCell);
ER_mean_vec   = cellfun(@mean,   ERRcell);

CT_median_vec = cellfun(@median, CTnegCell);
TW_median_vec = cellfun(@median, TWnegCell);
ER_median_vec = cellfun(@median, ERRcell);

pctFn         = @(v) prctile(v, bestPercentile);
CT_pX_vec     = cellfun(pctFn,   CTnegCell);
TW_pX_vec     = cellfun(pctFn,   TWnegCell);
ER_pX_vec     = cellfun(pctFn,   ERRcell);

% --- Best index per aggregation -----------------------------------------
[~, bestIdx_CT_mean]   = max(CT_mean_vec);
[~, bestIdx_TW_mean]   = max(TW_mean_vec);
[~, bestIdx_ER_mean]   = min(ER_mean_vec);

[~, bestIdx_CT_median] = max(CT_median_vec);
[~, bestIdx_TW_median] = max(TW_median_vec);
[~, bestIdx_ER_median] = min(ER_median_vec);

[~, bestIdx_CT_pX]     = max(CT_pX_vec);
[~, bestIdx_TW_pX]     = max(TW_pX_vec);
[~, bestIdx_ER_pX]     = min(ER_pX_vec);

% --- Select indices for figure highlighting (driven by BEST_AGG) --------
switch BEST_AGG
    case 'mean'
        bestIdx_CT = bestIdx_CT_mean;
        bestIdx_TW = bestIdx_TW_mean;
        bestIdx_ER = bestIdx_ER_mean;
    case 'median'
        bestIdx_CT = bestIdx_CT_median;
        bestIdx_TW = bestIdx_TW_median;
        bestIdx_ER = bestIdx_ER_median;
    case 'pX'
        bestIdx_CT = bestIdx_CT_pX;
        bestIdx_TW = bestIdx_TW_pX;
        bestIdx_ER = bestIdx_ER_pX;
    otherwise
        error('Unknown BEST_AGG: %s. Use ''mean'', ''median'', or ''pX''.', BEST_AGG);
end

% Specular codeword: phi1=0, phi2=0 — always index 1 with the shared codebook
specIdx = 1;

% Envelope min/max across all codewords
fMin_CT = min(cdf_CT);  fMax_CT = max(cdf_CT);
fMin_TW = min(cdf_TW);  fMax_TW = max(cdf_TW);
fMin_ER = min(cdf_ER);  fMax_ER = max(cdf_ER);

fprintf('Best codeword indices (%s, p=%d):  CT=%d  TW=%d  ER=%d\n\n', ...
        BEST_AGG, bestPercentile, bestIdx_CT, bestIdx_TW, bestIdx_ER);

%% ========================================================================
%% SECTION 6 — FIGURE GENERATION
%% ========================================================================
fprintf('=== Generating figures ===\n');

plotFn = @(xVals, cMat, cBase, fMin, fMax, labelStr, ...
           bestIdx, cFill, cBest, rawCell, rawDir, saveName) ...
    plotMetric_withLOS(xVals, cMat, cBase, fMin, fMax, labelStr, ...
                       bestIdx, specIdx, cFill, cBest, 'k', ...
                       rawCell, rawDir, saveName, figFolder, fileTag, nCodes);

plotFn(x_CT, cdf_CT, cdf_CT_dir, fMin_CT, fMax_CT, ...
       '$-\mathrm{Continuity}$',      bestIdx_CT, 'b', [0 0 1], CTnegCell, CTdir, 'NegCT');

plotFn(x_TW, cdf_TW, cdf_TW_dir, fMin_TW, fMax_TW, ...
       '$-\mathrm{Trustworthiness}$', bestIdx_TW, 'r', [1 0 0], TWnegCell, TWdir, 'NegTW');

plotFn(x_ER, cdf_ER, cdf_ER_dir, fMin_ER, fMax_ER, ...
       'Positioning Error (m)',        bestIdx_ER, 'g', [0 0.45 0], ERRcell, err_dir, 'PosErr');

%% ========================================================================
%% SECTION 7 — SAVE BEST-INDEX STRUCT  (lightweight append)
%%   All aggregation variants are stored so Main_DoubleEMS_Trajectory can
%%   load any combination without rerunning this script.
%%   codebook_all is also appended so the trajectory main can look up the
%%   phase values (phi1, phi2) for any chosen index.
%% ========================================================================
bestIdxStruct = struct( ...
    'CT_mean',    bestIdx_CT_mean,   'TW_mean',    bestIdx_TW_mean,   'ER_mean',    bestIdx_ER_mean,   ...
    'CT_median',  bestIdx_CT_median, 'TW_median',  bestIdx_TW_median, 'ER_median',  bestIdx_ER_median, ...
    'CT_pX',      bestIdx_CT_pX,     'TW_pX',      bestIdx_TW_pX,     'ER_pX',      bestIdx_ER_pX,     ...
    'percentile', bestPercentile,    'agg_used',   BEST_AGG);

save(dataFile, 'bestIdxStruct', 'codebook_all', '-append');

fprintf('Figures saved to:\n  %s\n\n', figFolder);
fprintf('=== Done ===\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LOCAL HELPER: ECDF figure with envelope + specular + best + baseline
%%
%%  Kept here (not in src/) because the specular/best highlighting and
%%  the mean-value text annotation are specific to this experiment.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function plotMetric_withLOS(xVals, cMat, cBase, fMin, fMax, labelStr, ...
                             bestIdx, specIdx, colorFill, colorBest, colorLOS, ...
                             rawCell, rawDir, saveName, figFolder, fileTag, nCodes)
% PLOTMETRIC_WITHLOS
% Draws one ECDF figure comparing the full static-codeword space against
% the direct-path-only baseline, and highlights the specular and best
% codewords.
%
% Inputs:
%   xVals      : [1 x nBins]   x-axis bin centres
%   cMat       : [nCodes x nBins]  CDF matrix, one row per codeword
%   cBase      : [1 x nBins]   CDF for the No-EMS baseline
%   fMin/fMax  : [1 x nBins]   point-wise min/max across codewords (envelope)
%   labelStr   : x-axis display label (may contain LaTeX)
%   bestIdx    : row index in cMat of the best codeword for this metric
%   specIdx    : row index of the specular (phi=0,0) codeword
%   colorFill  : shaded envelope colour (e.g. 'b', 'r', 'g')
%   colorBest  : RGB or named colour for the best-codeword line
%   colorLOS   : colour for the No-EMS baseline line
%   rawCell    : {nCodes x 1} cell of raw per-user metric vectors
%   rawDir     : vector of raw per-user metric values for the baseline
%   saveName   : clean filename prefix (no special characters)
%   figFolder  : output directory
%   fileTag    : file tag string appended to the saved filename
%   nCodes     : total number of codewords

    figure('Position', [100 100 800 600]);
    hold on; grid on; box on;

    % Shaded envelope (full codeword range)
    fill([xVals fliplr(xVals)], [fMin fliplr(fMax)], colorFill, ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
         'DisplayName', 'Codeword range');

    % All individual codewords as thin grey dashed (skip specular and best)
    for i = 1:nCodes
        if any(i == [specIdx, bestIdx]), continue; end
        plot(xVals, cMat(i,:), '--', 'LineWidth', 0.5, ...
             'Color', [0.6 0.6 0.6 0.6], 'HandleVisibility', 'off');
    end

    % Specular codeword (phi1=0, phi2=0)
    plot(xVals, cMat(specIdx,:), 'k-.', 'LineWidth', 2.0, ...
         'DisplayName', 'Specular ($\varphi\!=\!0,\,0$)');

    % Best codeword for this metric
    plot(xVals, cMat(bestIdx,:), '-', 'LineWidth', 2.5, ...
         'Color', colorBest, 'DisplayName', 'Best codeword');

    % No-EMS baseline
    plot(xVals, cBase, '-', 'LineWidth', 2.5, ...
         'Color', colorLOS, 'DisplayName', 'No EMS (direct only)');

    % Axis labels and formatting
    xlabel(labelStr, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('CDF',    'Interpreter', 'latex', 'FontSize', 18);
    title(['Double EMS: ' labelStr], 'Interpreter', 'latex', 'FontSize', 16);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'FontSize', 14);
    ylim([0 1]);
    if contains(labelStr, 'Error')
        xlim([0 150]);
    else
        xlim([min(xVals) 0]);
    end

    % Annotation: mean values of specular and best
    specAvg = mean(rawCell{specIdx});
    bestAvg = mean(rawCell{bestIdx});
    text(mean(xlim), 0.05, ...
         sprintf('\\textbf{Spec avg}=%.3f \\quad \\textbf{Best avg}=%.3f', ...
                 specAvg, bestAvg), ...
         'Interpreter', 'latex', 'FontSize', 14, ...
         'HorizontalAlignment', 'center');

    % Save
    savefig(fullfile(figFolder, [saveName '_' fileTag '.fig']));
    saveas(gcf,  fullfile(figFolder, [saveName '_' fileTag '.png']));
end

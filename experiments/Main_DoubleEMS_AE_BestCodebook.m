%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Main_DoubleEMS_AE_BestCodebook.m
%
%  Generates:
%    Fig. 7  — CDF of −Continuity and −Trustworthiness (AE, Double-EMS)
%              highlighting the specular codeword and the best codeword
%    Fig. 8  — CDF of Positioning Error (AE, Double-EMS)
%              highlighting the specular codeword and the best codeword
%
%  For each metric the figure shows:
%    • Shaded envelope   — range across all static codewords
%    • Thin grey dashed  — every individual codeword
%    • Specular (φ=0,0)  — reference codeword
%    • Best codeword     — selected by aggregation of choice
%    • No EMS baseline   — direct channel only
%
%  Embedding : semi-supervised Autoencoder (Python/TensorFlow)
%              Interface: MATLAB saves input_data.mat (R_Cov_BS, R_xy,
%              Supervision_percentage), Python writes ae_output.mat (xy_FCAE).
%
%  Same structure as Main_Double_EMS_BestCodebook.m.
%  No dissimilarity matrix is needed — the AE works directly on the
%  log-matrix features of the covariance matrices.
%
%  Run-or-load logic
%  -----------------
%  If a previously saved data file is found in outFolder the heavy
%  computation (codebook loop + AE calls + baseline) is skipped entirely
%  and the script jumps straight to the CDF computation and figure generation.
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

root  = get_project_root();
Data  = load_dataset(params);
state = prepare_data(Data, params);

%% ========================================================================
%% SECTION 1 — EXPERIMENT-SPECIFIC CHOICES  (override shared params here)
%% ========================================================================

% --- Supervision ---------------------------------------------------------
%     AE supervision is typically higher than t-SNE (more anchors needed)
Supervision_percentage = 0.15;
supervised_Len = round(params.Total_samples_used * Supervision_percentage);

% --- Best-codeword selection ---------------------------------------------
%     Same options as Main_Double_EMS_BestCodebook:
%       'mean' | 'median' | 'pX'  (bestPercentile only used with 'pX')
BEST_AGG       = 'pX';
bestPercentile = 90;

% --- Python environment and AE script ------------------------------------
%  PYTHON_MODE : 'venv'  → uses a local virtualenv folder (env_tensorflow_ae/)
%              : 'conda' → uses a conda environment named env_tensorflow_ae
%  Switch between the two by changing PYTHON_MODE only — nothing else changes.
PYTHON_MODE = 'venv';          % <-- change to 'conda' to use conda instead
ENV_NAME    = 'env_tensorflow_ae';   % name used for both venv folder and conda env

ae_script = fullfile(root, 'src', 'run_autoencoder.py');

% Build the shell command used to call Python (used in local_ae_eval).
% We cd to root first so that Python resolves input_data.mat and
% ae_output.mat relative to the project root (matching where MATLAB saves them).
if strcmp(PYTHON_MODE, 'venv')
    venv_python = fullfile(root, ENV_NAME, 'bin', 'python');
    ae_cmd = ['cd "' root '" && "' venv_python '" "' ae_script '"'];
else
    ae_cmd = ['cd "' root '" && conda run -n ' ENV_NAME ' python "' ae_script '"'];
end

% --- Temp files shared with Python (in project root) --------------------
ae_input_file  = fullfile(root, 'input_data.mat');
ae_output_file = fullfile(root, 'ae_output.mat');

% --- Codebook  (phi_options = deg2rad(0:4:40), 11 x 11 = 121 entries) ---
codebook     = generate_codebook(params);
nCodes       = codebook.nCodebooks;
codebook_all = codebook.all;          % [nCodes x 2]

% --- Derived convenience variables from state / params ------------------
Len     = state.Len;
R_xy    = state.R_xy;
ArrSize = params.Par.Tx.ArrSize;
N_real  = params.N_real;
M       = params.M;
N       = params.N;

%% ========================================================================
%% SECTION 2 — OUTPUT FOLDERS
%% ========================================================================

samplesStr = sprintf('%d_Samples', params.Total_samples_used);

figFolder = fullfile(root, 'Figures', samplesStr, 'DoubleEMS_AE', filesep);
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

fileTag = sprintf('DoubleEMS_AE_samp%d_super%.2f_nCodes%d', ...
                  params.Total_samples_used, Supervision_percentage, nCodes);

outFolder = fullfile(root, 'merged_results_doubleEMS_AE', fileTag, filesep);
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
    %% ENV SETUP — only runs once if the environment does not exist yet
    %% ----------------------------------------------------------------
    setup_script = fullfile(root, 'setup_env.sh');
    if strcmp(PYTHON_MODE, 'venv')
        env_exists = isfile(fullfile(root, ENV_NAME, 'bin', 'python'));
    else
        [env_check_status, ~] = system(['conda run -n ' ENV_NAME ' python --version']);
        env_exists = (env_check_status == 0);
    end

    if ~env_exists
        fprintf('=== Environment ''%s'' not found. Running setup_env.sh ===\n', ENV_NAME);
        status = system(['chmod +x "' setup_script '" && "' setup_script '" ' PYTHON_MODE ' ' ENV_NAME]);
        if status ~= 0
            error('Environment setup failed. Check the output above.');
        end
        fprintf('=== Environment ready. ===\n\n');
    else
        fprintf('=== Environment ''%s'' found. Skipping setup. ===\n\n', ENV_NAME);
    end

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

        % Build EMS-reflected channels
        H_ems1 = zeros(Len, ArrSize, size(state.H_i_EMS1, 3));
        H_ems2 = zeros(Len, ArrSize, size(state.H_i_EMS2, 3));
        for uu = 1:Len
            H_ems1(uu,:,:) = state.H_o_EMS1 * Theta1 * squeeze(state.H_i_EMS1(uu,:,:));
            H_ems2(uu,:,:) = state.H_o_EMS2 * Theta2 * squeeze(state.H_i_EMS2(uu,:,:));
        end
        H_tot = state.H_dir + H_ems1 + H_ems2;

        % Sample covariance matrices
        R_Cov_BS = local_compute_covariance(H_tot, Len, ArrSize, N_real);

        % AE embedding + CT/TW/ERR
        [CTneg, TWneg, err] = local_ae_eval( ...
            R_Cov_BS, R_xy, Supervision_percentage, supervised_Len, Len, ...
            ae_input_file, ae_output_file, ae_cmd, false);

        CTnegCell{idx} = CTneg;
        TWnegCell{idx} = TWneg;
        ERRcell{idx}   = err;
    end
    fprintf('Static codebook loop complete.\n\n');

    %% ----------------------------------------------------------------
    %% 3B — NO-EMS BASELINE  (direct channel only)
    %% ----------------------------------------------------------------
    fprintf('=== Computing No-EMS baseline ===\n');

    R_Cov_BS_dir = local_compute_covariance(state.H_dir, Len, ArrSize, N_real);
    [CTdir, TWdir, err_dir] = local_ae_eval( ...
        R_Cov_BS_dir, R_xy, Supervision_percentage, supervised_Len, Len, ...
        ae_input_file, ae_output_file, ae_cmd, true);

    fprintf('No-EMS baseline done.\n\n');

    %% ----------------------------------------------------------------
    %% 3C — SAVE
    %% ----------------------------------------------------------------
    fprintf('Saving data to:\n  %s\n', dataFile);
    save(dataFile, ...
        'CTnegCell', 'TWnegCell', 'ERRcell', ...
        'CTdir', 'TWdir', 'err_dir', ...
        'codebook_all', 'nCodes', 'supervised_Len', ...
        'params', 'Supervision_percentage', ...
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

tmp_CT = cellfun(@(v) histcounts(v, edgesCT, 'Norm','cdf'), CTnegCell, 'uni', 0);
tmp_TW = cellfun(@(v) histcounts(v, edgesCT, 'Norm','cdf'), TWnegCell, 'uni', 0);
tmp_ER = cellfun(@(v) histcounts(v, edgesER, 'Norm','cdf'), ERRcell,   'uni', 0);
cdf_CT = vertcat(tmp_CT{:});
cdf_TW = vertcat(tmp_TW{:});
cdf_ER = vertcat(tmp_ER{:});

cdf_CT_dir = histcounts(CTdir,    edgesCT, 'Norm','cdf');
cdf_TW_dir = histcounts(TWdir,    edgesCT, 'Norm','cdf');
cdf_ER_dir = histcounts(err_dir,  edgesER, 'Norm','cdf');

x_CT = edgesCT(2:end);
x_TW = x_CT;
x_ER = edgesER(2:end);

fprintf('CDFs computed.\n\n');

%% ========================================================================
%% SECTION 5 — BEST & SPECULAR CODEWORD INDICES
%%   All three aggregation methods computed; BEST_AGG selects which one is
%%   highlighted in figures.  CT/TW: negated → maximize.  ER: minimize.
%% ========================================================================

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

[~, bestIdx_CT_mean]   = max(CT_mean_vec);
[~, bestIdx_TW_mean]   = max(TW_mean_vec);
[~, bestIdx_ER_mean]   = min(ER_mean_vec);

[~, bestIdx_CT_median] = max(CT_median_vec);
[~, bestIdx_TW_median] = max(TW_median_vec);
[~, bestIdx_ER_median] = min(ER_median_vec);

[~, bestIdx_CT_pX]     = max(CT_pX_vec);
[~, bestIdx_TW_pX]     = max(TW_pX_vec);
[~, bestIdx_ER_pX]     = min(ER_pX_vec);

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

specIdx = 1;   % phi1=0, phi2=0 is always index 1

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
    plotMetric_withLOS_AE(xVals, cMat, cBase, fMin, fMax, labelStr, ...
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
%%   Stores all aggregation variants so a trajectory main can load any
%%   combination.  codebook_all appended for phi lookup.
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
%%  LOCAL FUNCTION: compute sample covariance matrices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function R_Cov_BS = local_compute_covariance(H_tot, Len, ArrSize, N_real)
% LOCAL_COMPUTE_COVARIANCE
%   Compute per-user sample covariance matrices from the channel H_tot.
%   R_Cov_BS : [Len x ArrSize x ArrSize]
    R_Cov_BS = zeros(Len, ArrSize, ArrSize);
    for uu = 1:Len
        H = squeeze(H_tot(uu, :, :));
        C = zeros(ArrSize);
        for r = 1:N_real
            Hs = H * (randn + 1j*randn) / sqrt(2);
            C  = C + Hs * Hs';
        end
        R_Cov_BS(uu, :, :) = C / N_real;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LOCAL FUNCTION: run AE and compute CT / TW / ERR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [CTneg, TWneg, err] = local_ae_eval( ...
        R_Cov_BS, R_xy, Supervision_percentage, supervised_Len, Len, ...
        ae_input_file, ae_output_file, ae_cmd, use_all_users)
% LOCAL_AE_EVAL
%   1. Saves R_Cov_BS + R_xy + Supervision_percentage to ae_input_file.
%   2. Runs the Python AE script via system(ae_cmd).
%   3. Loads xy_FCAE from ae_output_file.
%   4. Computes CT/TW and positioning error.
%
%   use_all_users = false : CT/TW/ERR computed on unlabeled users only
%                           (rows supervised_Len+1:end) — used for codebook loop
%   use_all_users = true  : CT/TW computed on ALL users (rows 1:end),
%                           ERR still on unlabeled only — used for No-EMS baseline
%                           (matches original code behaviour, "modified by reza")

    % --- Save inputs for Python ------------------------------------------
    save(ae_input_file, 'R_Cov_BS', 'R_xy', 'Supervision_percentage');

    % --- Run Python AE via shell -----------------------------------------
    status = system(ae_cmd);
    if status ~= 0
        error('AE Python script failed (exit code %d). Check terminal output.', status);
    end

    % --- Load AE output --------------------------------------------------
    ae_out  = load(ae_output_file, 'xy_FCAE');
    xy_FCAE = ae_out.xy_FCAE;   % [Len x 2], already in physical units

    % --- Positioning error (unlabeled users only, always) ----------------
    err = sqrt(sum((xy_FCAE - R_xy).^2, 2));
    err = err(supervised_Len+1:end);

    % --- CT / TW ---------------------------------------------------------
    KK    = 0.05 * Len;
    DistT = squareform(pdist(R_xy));
    DistZ = squareform(pdist(xy_FCAE));
    TC    = drscore_pointwise(DistT, DistZ, KK, 'TC');

    if use_all_users
        % No-EMS baseline: use all users (matches original "modified by reza")
        CTneg = -TC(1:end, 1);
        TWneg = -TC(1:end, 2);
    else
        % Codebook loop: unlabeled users only
        CTneg = -TC(supervised_Len+1:end, 1);
        TWneg = -TC(supervised_Len+1:end, 2);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LOCAL HELPER: ECDF figure with envelope + specular + best + baseline
%%
%%  Identical to plotMetric_withLOS in Main_Double_EMS_BestCodebook except
%%  the title says "Double EMS with AE".
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotMetric_withLOS_AE(xVals, cMat, cBase, fMin, fMax, labelStr, ...
                                 bestIdx, specIdx, colorFill, colorBest, colorLOS, ...
                                 rawCell, rawDir, saveName, figFolder, fileTag, nCodes)

    figure('Position', [100 100 800 600]);
    hold on; grid on; box on;

    % Shaded envelope
    fill([xVals fliplr(xVals)], [fMin fliplr(fMax)], colorFill, ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
         'DisplayName', 'Codeword range');

    % All individual codewords (thin grey dashed, skip specular and best)
    for i = 1:nCodes
        if any(i == [specIdx, bestIdx]), continue; end
        plot(xVals, cMat(i,:), '--', 'LineWidth', 0.5, ...
             'Color', [0.6 0.6 0.6 0.6], 'HandleVisibility', 'off');
    end

    % Specular codeword
    plot(xVals, cMat(specIdx,:), 'k-.', 'LineWidth', 2.0, ...
         'DisplayName', 'Specular ($\varphi\!=\!0,\,0$)');

    % Best codeword
    plot(xVals, cMat(bestIdx,:), '-', 'LineWidth', 2.5, ...
         'Color', colorBest, 'DisplayName', 'Best codeword');

    % No-EMS baseline
    plot(xVals, cBase, '-', 'LineWidth', 2.5, ...
         'Color', colorLOS, 'DisplayName', 'No EMS (direct only)');

    xlabel(labelStr, 'Interpreter', 'latex', 'FontSize', 18);
    ylabel('CDF',    'Interpreter', 'latex', 'FontSize', 18);
    title(['Double EMS with AE: ' labelStr], 'Interpreter', 'latex', 'FontSize', 16);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);
    set(gca, 'FontSize', 14);
    ylim([0 1]);
    if contains(labelStr, 'Error')
        xlim([0 150]);
    else
        xlim([min(xVals) 0]);
    end

    specAvg = mean(rawCell{specIdx});
    bestAvg = mean(rawCell{bestIdx});
    text(mean(xlim), 0.05, ...
         sprintf('\\textbf{Spec avg}=%.3f \\quad \\textbf{Best avg}=%.3f', ...
                 specAvg, bestAvg), ...
         'Interpreter', 'latex', 'FontSize', 14, ...
         'HorizontalAlignment', 'center');

    savefig(fullfile(figFolder, [saveName '_' fileTag '.fig']));
    saveas(gcf,  fullfile(figFolder, [saveName '_' fileTag '.png']));
end

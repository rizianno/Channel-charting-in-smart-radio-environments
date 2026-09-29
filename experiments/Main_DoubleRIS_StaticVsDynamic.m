%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Main_DoubleEMS_StaticVsDynamic.m
%
%  Generates:
%    Fig. 7  — CDF of −Continuity and −Trustworthiness (t-SNE, Double-EMS)
%    Fig. 8  — CDF of Positioning Error in metres        (t-SNE, Double-EMS)
%    Fig. 10 — CDF of SNR and CMD Dissimilarity
%
%  Compares across:
%    • Static codebook EMSs  — full envelope (grey fill) + individual curves
%    • No EMS / LOS baseline
%    • Dynamic   (|H_o| * |H_i| incoherent combination)
%    • Dynamic |H_dir|  (same but with |H_dir| on direct path)
%    • Random Θ  (independent random phase per element, per user)
%
%  Dissimilarity : CMD  (Correlation Matrix Distance)
%  Embedding     : semi-supervised t-SNE
%
%  Run-or-load logic
%  -----------------
%  If a previously saved data file is found in outFolder the heavy
%  computation is skipped entirely and the script jumps straight to the
%  CDF computation and figure generation.  This lets you re-plot or
%  adjust figure aesthetics without rerunning the simulation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% ------------------------------------------------------------------------
% INIT
%% ------------------------------------------------------------------------

if isempty(mfilename('fullpath'))
    % Running from command window → fallback
    root = pwd;
else
    root = fileparts(mfilename('fullpath'));
end

addpath(genpath(fullfile(root, '..', 'src')));
params = init_parameters();

%% ------------------------------------------------------------------------
% LOAD DATA
%% ------------------------------------------------------------------------
root = get_project_root();
resultsFolder = fullfile(root, 'results');

Data  = load_dataset(params);
state = prepare_data(Data, params);

%% ========================================================================
%% SECTION 1 — USER CHOICES  (edit here)
%% ========================================================================

% --- Dataset (update path for your machine) ---
dataPath = fullfile('C:', 'Users', 'Reza Aghazadeh', ...
    'OneDrive - Politecnico di Milano', 'RIS_CHannel_for_CC', ...
    'HRIS_BS30janLOSNLOS__uniform32_8_RISy4z5751bounceuser_3200_23tx_power.mat');

% --- Experiment settings ---
dissim_type            = 'LE';
Total_samples_used     = 3200;
Supervision_percentage = 0.30;
supervised_Len         = round(Total_samples_used * Supervision_percentage);

% --- EMS / codebook ---
M = 60;  N = 60;
phi_options_1 = deg2rad(0:4:4);   % phase gradient options for EMS 1
phi_options_2 = deg2rad(0:4:4);   % phase gradient options for EMS 2
nCodes        = numel(phi_options_1) * numel(phi_options_2);
codebook_all  = combvec(phi_options_1, phi_options_2).';  % [nCodes x 2]

% --- System ---
ArrSize = 32;      % BS antenna count
N_real  = 500;     % noise realisations for sample covariance
snr_tx  = 25 - (-90);   % TX SNR offset [dB]

%% ========================================================================
%% SECTION 2 — OUTPUT FOLDERS
%% ========================================================================

% % % % projectRoot = fullfile(thisDir, '..');
samplesStr  = sprintf('%d_Samples', Total_samples_used);

figFolder = fullfile(root, 'Figures', samplesStr, ...
            ['DoubleEMS_ScaledV1_' dissim_type], filesep);
if ~exist(figFolder, 'dir'), mkdir(figFolder); end

fileTag = sprintf('DoubleEMS_ScaledV1_%s_samp%d_super%.2f_nCodes%d', ...
                  dissim_type, Total_samples_used, ...
                  Supervision_percentage, nCodes);

outFolder = fullfile(root, 'merged_results_doubleEMS_ScaledV1', ...
                     fileTag, filesep);
if ~exist(outFolder, 'dir'), mkdir(outFolder); end

%% ========================================================================
%% SECTION 3 — RUN OR LOAD
%% ========================================================================

dataFile = fullfile(outFolder, ['MergedData_' fileTag '.mat']);

if isfile(dataFile)
    % ------------------------------------------------------------------
    %  Fast path: load previously saved results and skip all computation
    % ------------------------------------------------------------------
    fprintf('=== Saved data found ===\n');
    fprintf('Loading from:\n  %s\n\n', dataFile);
    load(dataFile);
    Len = Total_samples_used;
    fprintf('Load complete. Skipping computation.\n\n');

else
    % ------------------------------------------------------------------
    %  Full computation path
    % ------------------------------------------------------------------
    fprintf('=== No saved data found. Running full computation. ===\n\n');

    %% ----------------------------------------------------------------
    %% 3A — LOAD DATASET
    %% ----------------------------------------------------------------
    fprintf('Loading dataset ...\n');
    % % % % % % % assert(isfile(dataPath), ...
    % % % % % % %     ['Dataset not found:\n  %s\n' ...
    % % % % % % %      'Update the dataPath variable in SECTION 1.'], dataPath);
    % % % % % % % 
    % % % % % % % Data = load(dataPath);

    Hi_1 = Data.H_I_RIS1;   Ho_1 = Data.H_o_RIS1;
    Hi_2 = Data.H_I_RIS2;   Ho_2 = Data.H_o_RIS2;

    idxUsers = 1:Total_samples_used;
    H_dir    = Data.H_LOS_1bounce(idxUsers, :, :);
    R_xy     = Data.R_xy(idxUsers, :);
    Len      = numel(idxUsers);
    Hi_1     = Hi_1(idxUsers, :, :);
    Hi_2     = Hi_2(idxUsers, :, :);

    fprintf('Dataset loaded: %d users.\n\n', Len);

    %% ----------------------------------------------------------------
    %% 3B — STATIC CODEBOOK LOOP
    %% ----------------------------------------------------------------
    fprintf('=== Static codebook loop (%d entries) ===\n', nCodes);

    CTnegCell = cell(nCodes, 1);   TWnegCell = cell(nCodes, 1);
    ERRcell   = cell(nCodes, 1);   GainCell  = cell(nCodes, 1);
    DvecCell  = cell(nCodes, 1);   SNRCell   = cell(nCodes, 1);

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

        % Build static EMS-reflected channels
        H_ris1 = zeros(Len, ArrSize, size(Hi_1, 3));
        H_ris2 = zeros(Len, ArrSize, size(Hi_2, 3));
        for uu = 1:Len
            H_ris1(uu,:,:) = Ho_1 * Theta1 * squeeze(Hi_1(uu,:,:));
            H_ris2(uu,:,:) = Ho_2 * Theta2 * squeeze(Hi_2(uu,:,:));
        end
        H_tot = sqrt(8) * (H_dir + 4*H_ris1 + 4*H_ris2);

        % Evaluate charting metrics via src utility
        [CTn, TWn, ERR, gainVec, dvec_dB] = evaluate_channel_configuration( ...
            H_tot, R_xy, supervised_Len, ArrSize, N_real, ...
            @CompDismilarity_CMD);

        CTnegCell{idx} = CTn;
        TWnegCell{idx} = TWn;
        ERRcell{idx}   = ERR;
        GainCell{idx}  = gainVec;
        DvecCell{idx}  = dvec_dB;
        SNRCell{idx}   = gainVec + snr_tx;
    end
    fprintf('Static codebook loop complete.\n\n');

    %% ----------------------------------------------------------------
    %% 3C — BASELINES  (dynamic / random / direct-only)
    %% ----------------------------------------------------------------
    fprintf('=== Computing baselines ===\n');

    % Shared evaluation wrapper (CMD dissimilarity, same settings)
    evalFn = @(Htot) evaluate_channel_configuration( ...
                Htot, R_xy, supervised_Len, ArrSize, N_real, ...
                @CompDismilarity_CMD);

    % --- Precompute incoherent EMS channels (reused by dyn #1 and #2) ---
    H_r1_dyn = zeros(Len, ArrSize, size(Hi_1, 3));
    H_r2_dyn = zeros(Len, ArrSize, size(Hi_2, 3));
    for uu = 1:Len
        H_r1_dyn(uu,:,:) = abs(Ho_1) * abs(squeeze(Hi_1(uu,:,:)));
        H_r2_dyn(uu,:,:) = abs(Ho_2) * abs(squeeze(Hi_2(uu,:,:)));
    end

    % Dynamic #1 — incoherent (|H_o| * |H_i|) combination
    H_dyn = sqrt(8) * (H_dir + 4*H_r1_dyn + 4*H_r2_dyn);
    [CT_dyn, TW_dyn, ERR_dyn, G_dyn, D_dyn] = evalFn(H_dyn);
    SNR_dyn = G_dyn + snr_tx;
    fprintf('  Dynamic baseline done.\n');

    % Dynamic #2 — same but also take |H_dir|
    H_dynAbs = sqrt(8) * (abs(H_dir) + 4*H_r1_dyn + 4*H_r2_dyn);
    [CT_dyn2, TW_dyn2, ERR_dyn2, G_dyn2, D_dyn2] = evalFn(H_dynAbs);
    SNR_dyn2 = G_dyn2 + snr_tx;
    fprintf('  Dynamic |H_dir| baseline done.\n');

    % Random phases — independent uniform phase per element, per user
    rng('default');
    H_r1_rnd = zeros(Len, ArrSize, size(Hi_1, 3));
    H_r2_rnd = zeros(Len, ArrSize, size(Hi_2, 3));
    for uu = 1:Len
        Theta1r = diag(exp(-1j * 2*pi * rand(M*N, 1)));
        Theta2r = diag(exp(-1j * 2*pi * rand(M*N, 1)));
        H_r1_rnd(uu,:,:) = Ho_1 * Theta1r * squeeze(Hi_1(uu,:,:));
        H_r2_rnd(uu,:,:) = Ho_2 * Theta2r * squeeze(Hi_2(uu,:,:));
    end
    H_rnd = sqrt(8) * (H_dir + 4*H_r1_rnd + 4*H_r2_rnd);
    [CT_rnd, TW_rnd, ERR_rnd, G_rnd, D_rnd] = evalFn(H_rnd);
    SNR_rnd = G_rnd + snr_tx;
    fprintf('  Random-phase baseline done.\n');

    % Direct-path-only (no EMS / LOS) baseline
    H_dir_scaled = sqrt(8) * H_dir;
    [CTdir, TWdir, err_dir, gain_dir, dvec_dir] = evalFn(H_dir_scaled);
    SNR_dir = gain_dir + snr_tx;
    fprintf('  Direct (LOS) baseline done.\n\n');

    %% ----------------------------------------------------------------
    %% 3D — SAVE  (all raw per-user metric vectors)
    %% ----------------------------------------------------------------
    fprintf('Saving data to:\n  %s\n', dataFile);
    save(dataFile, ...
        'CTnegCell', 'TWnegCell', 'ERRcell', ...
        'GainCell',  'DvecCell',  'SNRCell', ...
        'CT_dyn',  'TW_dyn',  'ERR_dyn',  'G_dyn',  'D_dyn',  'SNR_dyn',  ...
        'CT_dyn2', 'TW_dyn2', 'ERR_dyn2', 'G_dyn2', 'D_dyn2', 'SNR_dyn2', ...
        'CT_rnd',  'TW_rnd',  'ERR_rnd',  'G_rnd',  'D_rnd',  'SNR_rnd',  ...
        'CTdir',   'TWdir',   'err_dir',  'gain_dir','dvec_dir','SNR_dir', ...
        'R_xy', 'codebook_all', 'phi_options_1', 'phi_options_2', ...
        'snr_tx', 'Total_samples_used', 'Supervision_percentage', ...
        'supervised_Len', '-v7.3');
    fprintf('Data saved successfully.\n\n');
end

%% ========================================================================
%% SECTION 4 — EMPIRICAL CDF COMPUTATION
%%   (always runs — fast, so no need to cache; re-bin on every launch)
%% ========================================================================
fprintf('=== Computing empirical CDFs ===\n');

edgesCT = linspace(-1, 0, 201);
edgesER = linspace(0, 150, 201);

gain_all_dB = [cell2mat(GainCell(:)) ; gain_dir(:) ; ...
               G_dyn(:) ; G_dyn2(:) ; G_rnd(:)];
edgesG   = linspace(min(gain_all_dB), max(gain_all_dB), 201);
edges_snr = edgesG + snr_tx;

edgesD = linspace(-140, 0, 201);

x_CT  = edgesCT(2:end);
x_TW  = x_CT;
x_ER  = edgesER(2:end);
x_snr = edges_snr(2:end);
x_D   = edgesD(2:end);

% --- Static codebook CDFs (cell arrays, one entry per codeword) ----------
cdf_CT  = cellfun(@(v) histcounts(v, edgesCT,   'Norm','cdf'), CTnegCell,'uni',0);
cdf_TW  = cellfun(@(v) histcounts(v, edgesCT,   'Norm','cdf'), TWnegCell,'uni',0);
cdf_ER  = cellfun(@(v) histcounts(v, edgesER,   'Norm','cdf'), ERRcell,  'uni',0);
cdf_SNR = cellfun(@(v) histcounts(v, edges_snr, 'Norm','cdf'), SNRCell,  'uni',0);
cdf_D   = cellfun(@(v) histcounts(v, edgesD,    'Norm','cdf'), DvecCell, 'uni',0);

% --- Direct (LOS) baseline -----------------------------------------------
cdf_CT_dir  = histcounts(CTdir,    edgesCT,   'Norm','cdf');
cdf_TW_dir  = histcounts(TWdir,    edgesCT,   'Norm','cdf');
cdf_ER_dir  = histcounts(err_dir,  edgesER,   'Norm','cdf');
cdf_snr_dir = histcounts(SNR_dir,  edges_snr, 'Norm','cdf');
cdf_D_dir   = histcounts(dvec_dir, edgesD,    'Norm','cdf');

% --- Dynamic baseline #1 -------------------------------------------------
cdf_CT_dyn  = histcounts(CT_dyn,  edgesCT,   'Norm','cdf');
cdf_TW_dyn  = histcounts(TW_dyn,  edgesCT,   'Norm','cdf');
cdf_ER_dyn  = histcounts(ERR_dyn, edgesER,   'Norm','cdf');
cdf_snr_dyn = histcounts(SNR_dyn, edges_snr, 'Norm','cdf');
cdf_D_dyn   = histcounts(D_dyn,   edgesD,    'Norm','cdf');

% --- Dynamic baseline #2 (|H_dir|) ---------------------------------------
cdf_CT_dyn2  = histcounts(CT_dyn2,  edgesCT,   'Norm','cdf');
cdf_TW_dyn2  = histcounts(TW_dyn2,  edgesCT,   'Norm','cdf');
cdf_ER_dyn2  = histcounts(ERR_dyn2, edgesER,   'Norm','cdf');
cdf_snr_dyn2 = histcounts(SNR_dyn2, edges_snr, 'Norm','cdf');
cdf_D_dyn2   = histcounts(D_dyn2,   edgesD,    'Norm','cdf');

% --- Random-phase baseline -----------------------------------------------
cdf_CT_rnd  = histcounts(CT_rnd,  edgesCT,   'Norm','cdf');
cdf_TW_rnd  = histcounts(TW_rnd,  edgesCT,   'Norm','cdf');
cdf_ER_rnd  = histcounts(ERR_rnd, edgesER,   'Norm','cdf');
cdf_snr_rnd = histcounts(SNR_rnd, edges_snr, 'Norm','cdf');
cdf_D_rnd   = histcounts(D_rnd,   edgesD,    'Norm','cdf');

fprintf('CDFs computed.\n\n');

%% ========================================================================
%% SECTION 5 — FIGURE GENERATION
%% ========================================================================
fprintf('=== Generating figures ===\n');

% Convenience wrapper:
%   plotFn(x, displayLabel, fileLabel, staticCDFs, cdfsLOS, cdfs_dyn, ...)
plotFn = @(x, lab, saveName, cStat, cLOS, cDyn, cDyn2, cRnd, cFill) ...
    local_plotECDF(x, lab, saveName, cStat, cLOS, cDyn, cDyn2, cRnd, ...
                   cFill, figFolder, fileTag);

plotFn(x_CT,      '$-\mathrm{Continuity}$',      'NegCT', ...
       cdf_CT,  cdf_CT_dir,  cdf_CT_dyn,  cdf_CT_dyn2,  cdf_CT_rnd,  'b');

plotFn(x_TW,      '$-\mathrm{Trustworthiness}$', 'NegTW', ...
       cdf_TW,  cdf_TW_dir,  cdf_TW_dyn,  cdf_TW_dyn2,  cdf_TW_rnd,  'r');

plotFn(x_ER,      'Positioning Error (m)',        'PosErr', ...
       cdf_ER,  cdf_ER_dir,  cdf_ER_dyn,  cdf_ER_dyn2,  cdf_ER_rnd,  'g');

% NOTE: x_snr+10 preserves the manual axis offset present in the
%       original monolithic script (shifts x-axis by 10 dB for display).
plotFn(x_snr+10,  'SNR (dB)',                     'SNR', ...
       cdf_SNR, cdf_snr_dir, cdf_snr_dyn, cdf_snr_dyn2, cdf_snr_rnd, 'm');

plotFn(x_D,       'CMD Dissimilarity (dB)',        'Dissim', ...
       cdf_D,   cdf_D_dir,   cdf_D_dyn,   cdf_D_dyn2,   cdf_D_rnd,   'c');

fprintf('Figures saved to:\n  %s\n\n', figFolder);
fprintf('=== Done ===\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  LOCAL HELPER: ECDF figure with static envelope + 4 named curves
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function local_plotECDF(x, lab, saveName, cdfStat, ...
                         cLOS, cDyn, cDyn2, cRnd, ...
                         cFill, figFolder, fTag)
% LOCAL_PLOTECDF
% Draws a single ECDF figure with:
%   - a shaded envelope for the static codebook range
%   - thin grey dashed lines for every individual codeword
%   - four named solid curves: LOS, Dynamic, Dynamic|H_dir|, Random Θ
%
% Inputs:
%   x        : [1 x 200] x-axis bin centres
%   lab      : display label for x-axis (may include LaTeX)
%   saveName : clean string used in the saved file name (no spaces / $)
%   cdfStat  : {nCodes x 1} cell array of CDF row vectors (static codewords)
%   cLOS     : [1 x 200] CDF — No EMS / LOS baseline
%   cDyn     : [1 x 200] CDF — Dynamic (incoherent)
%   cDyn2    : [1 x 200] CDF — Dynamic |H_dir|
%   cRnd     : [1 x 200] CDF — Random phases
%   cFill    : colour string for the static envelope fill (e.g. 'b', 'r')
%   figFolder: path where figures are saved
%   fTag     : file tag string appended to the figure name

    cdfMat = vertcat(cdfStat{:});          % [nCodes x nBins]
    fMin   = min(cdfMat, [], 1);
    fMax   = max(cdfMat, [], 1);

    figure('Position', [90 90 820 600]);
    hold on; box on; grid on;

    % Shaded static envelope
    fill([x fliplr(x)], [fMin fliplr(fMax)], cFill, ...
         'FaceAlpha', 0.15, 'EdgeColor', 'none', ...
         'DisplayName', 'Static codeword range');

    % Individual static codewords (grey dashed, excluded from legend)
    for k = 1:size(cdfMat, 1)
        plot(x, cdfMat(k,:), '--', 'Color', [0.7 0.7 0.7], ...
             'HandleVisibility', 'off');
    end

    % Named curves
    plot(x, cLOS,  'k-',  'LineWidth', 2.4, 'DisplayName', 'No EMS (LOS)');
    plot(x, cDyn,  'b-',  'LineWidth', 2.0, 'DisplayName', 'Dynamic');
    plot(x, cDyn2, 'g-',  'LineWidth', 2.0, 'DisplayName', 'Dynamic $|H_{\mathrm{dir}}|$');
    plot(x, cRnd,  'm-',  'LineWidth', 2.0, 'DisplayName', 'Random $\Theta$');

    xlabel(lab,   'Interpreter', 'latex', 'FontSize', 18);
    ylabel('CDF', 'Interpreter', 'latex', 'FontSize', 18);
    title(['Double EMS — ' lab], 'Interpreter', 'latex', 'FontSize', 16);
    legend('Location', 'best', 'Interpreter', 'latex', 'FontSize', 14);

    ylim([0 1]);
    set(gca, 'FontSize', 14);

    % Sensible x-limits per metric type
    if contains(lab, 'Error')
        xlim([0 150]);
    elseif contains(lab, '-')       % negated metrics live in [-1, 0]
        xlim([min(x) 0]);
    else
        xlim([min(x) max(x)]);
    end

    % Save as .fig using the clean saveName (avoids special chars)
    savefig(fullfile(figFolder, [saveName '_' fTag '.fig']));

end

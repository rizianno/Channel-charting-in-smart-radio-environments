%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  DOUBLE-RIS  (Windows paths, −CT/−TW, 90-pct minimisation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;

%% PYTHON / ENVIRONMENT CONFIG  -------------------------------------------
PYTHON_MODE = 'venv';           % 'venv' | 'conda'
ENV_NAME    = 'env_tensorflow_ae';

%% PROJECT ROOT  ----------------------------------------------------------
% This file lives in  <root>/src/  so go up one level to get root
src_dir    = fileparts(mfilename('fullpath'));
root       = fileparts(src_dir);

% Paths used for the Python interface
ae_script      = fullfile(src_dir, 'run_autoencoder.py');
ae_input_file  = fullfile(root, 'input_data.mat');
ae_output_file = fullfile(root, 'ae_output.mat');

% Build the shell command that activates the env and runs the AE script
if strcmp(PYTHON_MODE, 'venv')
    venv_python = fullfile(root, ENV_NAME, 'bin', 'python');
    ae_cmd = ['cd "' root '" && "' venv_python '" "' ae_script '"'];
else
    ae_cmd = ['cd "' root '" && conda run -n ' ENV_NAME ' python "' ae_script '"'];
end

% Auto-setup Python environment if it does not exist yet
setup_script = fullfile(root, 'setup_env.sh');
if strcmp(PYTHON_MODE, 'venv')
    env_exists = isfile(fullfile(root, ENV_NAME, 'bin', 'python'));
else
    [env_check_status, ~] = system(['conda run -n ' ENV_NAME ' python --version']);
    env_exists = (env_check_status == 0);
end
if ~env_exists
    fprintf('Python environment not found. Running setup_env.sh ...\n');
    status = system(['chmod +x "' setup_script '" && "' setup_script '" ' PYTHON_MODE ' ' ENV_NAME]);
    if status ~= 0
        error('setup_env.sh failed (exit code %d). Check the script.', status);
    end
end

% Make sure src functions (drscore_pointwise, etc.) are on the MATLAB path
addpath(src_dir);

%% USER CHOICES -----------------------------------------------------------
dissim_type            = 'LE';      % 'LE' or 'CMD'
Total_samples_used     = 3200;
Supervision_percentage = 0.30;
supervised_Len         = round(Total_samples_used*Supervision_percentage);

%% LOAD DATASET -----------------------------------------------------------
dataPath = fullfile(root, 'data', 'HEMS_dataset.mat');
assert(isfile(dataPath), 'Dataset not found at: %s\nCheck README.', dataPath);
Data  = load(dataPath);
Hi_1  = Data.H_I_RIS1;  Ho_1 = Data.H_o_RIS1;
Hi_2  = Data.H_I_RIS2;  Ho_2 = Data.H_o_RIS2;

%% RANDOM USER SUBSET -----------------------------------------------------
idxUsers   = 1:Total_samples_used;                    % or  randi(...)
H_dir      = Data.H_LOS_1bounce(idxUsers,:,:);
R_xy       = Data.R_xy(idxUsers,:);
Len        = length(idxUsers);
Hi_1       = Hi_1(idxUsers,:,:);   Hi_2 = Hi_2(idxUsers,:,:);

%% CONSTANTS --------------------------------------------------------------
M = 60;  N = 60;
% % % % % phi_options_1 = deg2rad(0:4:40);
% % % % % phi_options_2 = deg2rad(0:4:40);

phi_options_1 = deg2rad(0:4:4);
phi_options_2 = deg2rad(0:4:4);

nCodes   = numel(phi_options_1)*numel(phi_options_2);
codebook_all = combvec(phi_options_1,phi_options_2).';   % each row = [φ1 φ2]

%% OUTPUT FOLDERS ---------------------------------------------------------
samplesStr = sprintf('%d_Samples',Total_samples_used);
figFolder  = [fullfile(root,'Figures',samplesStr,['AE_DoubleRIS_' dissim_type]) filesep];
if ~exist(figFolder,'dir'), mkdir(figFolder); end
fileTag    = sprintf('DoubleRIS_%s_samp%d_super%.2f_nCodes%d', ...
                     dissim_type,Total_samples_used,Supervision_percentage,nCodes);
outFolder  = [fullfile(root,['Results_AE_doubleRIS_' fileTag]) filesep];
if ~exist(outFolder,'dir'), mkdir(outFolder); end

%% PRE-ALLOC --------------------------------------------------------------
CTnegCell = cell(nCodes,1);   TWnegCell = cell(nCodes,1);   ERRcell = cell(nCodes,1);

%% MAIN LOOP --------------------------------------------------------------
Par.Tx.ArrSize = 32;    N_real = 500;
for idx = 1:nCodes
    phi1 = codebook_all(idx,1);   phi2 = codebook_all(idx,2);
    fprintf('[%3d/%3d]  φ1=%.2f  φ2=%.2f  (rad)\n',idx,nCodes,phi1,phi2);

    % -------- Channel construction (identical to Linux script) ----------
    phase1 = kron(phi1*(1:M),ones(1,N));
    phase2 = kron(phi2*(1:M),ones(1,N));
    Theta1 = diag(exp(-1j*phase1));   Theta2 = diag(exp(-1j*phase2));

    H_ris1 = zeros(Len,Par.Tx.ArrSize,size(Hi_1,3));
    H_ris2 = H_ris1;
    for uu = 1:Len
        H_ris1(uu,:,:) = Ho_1 * Theta1 * squeeze(Hi_1(uu,:,:));
        H_ris2(uu,:,:) = Ho_2 * Theta2 * squeeze(Hi_2(uu,:,:));
    end
    H_tot = H_dir + H_ris1 + H_ris2;

    % -------- Covariance, dissimilarity, t-SNE, metrics -----------------
    R_Cov_BS = zeros(Len,Par.Tx.ArrSize,Par.Tx.ArrSize);
    for uu = 1:Len
        H = squeeze(H_tot(uu,:,:));
        C = zeros(Par.Tx.ArrSize);
        for r = 1:N_real
            Hs = H*(randn+1j*randn)/sqrt(2);
            C  = C + Hs*Hs';
        end
        R_Cov_BS(uu,:,:) = C/N_real;
    end

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   % Save the data to a .mat file (absolute path; include Supervision_percentage for Python)
   save(ae_input_file, 'R_Cov_BS', 'R_xy', 'Supervision_percentage');

   % Call the Python script via the configured environment
   status = system(ae_cmd);
   if status ~= 0
       error('AE Python script failed (exit code %d).', status);
   end

   % Load the output from Python
   data = load(ae_output_file);
   xy_FCAE = data.xy_FCAE;

   disp('Latent output from Python:')
   % disp(latent)
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    err  = sqrt(sum((xy_FCAE-R_xy).^2,2));
    err  = err(supervised_Len+1:end);

    KK   = 0.05*Len;
    DistT  = squareform(pdist(R_xy));
    DistZ  = squareform(pdist(xy_FCAE));
    TC     = drscore_pointwise(DistT,DistZ,KK,'TC');
    CTneg  = -TC(supervised_Len+1:end,1);
    TWneg  = -TC(supervised_Len+1:end,2);

    CTnegCell{idx}=CTneg;   TWnegCell{idx}=TWneg;   ERRcell{idx}=err;
end

%% NO-EMS BASELINE --------------------------------------------------------
fprintf('=== Computing No-EMS baseline (direct channel only) ===\n');

R_Cov_BS_dir = zeros(Len, Par.Tx.ArrSize, Par.Tx.ArrSize);
for uu = 1:Len
    H = squeeze(H_dir(uu,:,:));
    C = zeros(Par.Tx.ArrSize);
    for r = 1:N_real
        Hs = H*(randn+1j*randn)/sqrt(2);
        C  = C + Hs*Hs';
    end
    R_Cov_BS_dir(uu,:,:) = C/N_real;
end

R_Cov_BS = R_Cov_BS_dir;   % reuse variable name expected by Python script
save(ae_input_file, 'R_Cov_BS', 'R_xy', 'Supervision_percentage');
status = system(ae_cmd);
if status ~= 0
    error('AE Python script (No-EMS baseline) failed (exit code %d).', status);
end
data_dir    = load(ae_output_file);
xy_FCAE_dir = data_dir.xy_FCAE;

err_dir  = sqrt(sum((xy_FCAE_dir - R_xy).^2, 2));
err_dir  = err_dir(supervised_Len+1:end);

KK_dir   = 0.05*Len;
DistT_dir = squareform(pdist(R_xy));
DistZ_dir = squareform(pdist(xy_FCAE_dir));
TC_dir    = drscore_pointwise(DistT_dir, DistZ_dir, KK_dir, 'TC');
CTdir     = -TC_dir(supervised_Len+1:end, 1);
TWdir     = -TC_dir(supervised_Len+1:end, 2);

fprintf('No-EMS baseline done.\n\n');

%% BUILD CDF MATRICES (identical logic) -----------------------------------
edgesCT = linspace(-1,0,201);   edgesER = linspace(0,150,201);
cdf_CT  = cellfun(@(x)histcounts(x,edgesCT,'norm','cdf'),CTnegCell,'uni',0);
cdf_TW  = cellfun(@(x)histcounts(x,edgesCT,'norm','cdf'),TWnegCell,'uni',0);
cdf_ER  = cellfun(@(x)histcounts(x,edgesER,'norm','cdf'),ERRcell ,'uni',0);
cdf_CT  = vertcat(cdf_CT{:});   cdf_TW = vertcat(cdf_TW{:});   cdf_ER = vertcat(cdf_ER{:});
x_CT = edgesCT(2:end);  x_ER = edgesER(2:end);  x_TW = x_CT;

% Baseline CDFs (single row vectors)
cdf_CT_dir = histcounts(CTdir,   edgesCT, 'Norm', 'cdf');
cdf_TW_dir = histcounts(TWdir,   edgesCT, 'Norm', 'cdf');
cdf_ER_dir = histcounts(err_dir, edgesER, 'Norm', 'cdf');

%% "BEST" & "SPECULAR" ----------------------------------------------------
pctVal = @(cdf,x,p) x(find(cdf>=p,1,'first'));
CT90 = arrayfun(@(i)pctVal(cdf_CT(i,:),x_CT,0.9),1:nCodes);
TW90 = arrayfun(@(i)pctVal(cdf_TW(i,:),x_TW,0.9),1:nCodes);
ER90 = arrayfun(@(i)pctVal(cdf_ER(i,:),x_ER,0.9),1:nCodes);

[~,bestIdx_CT] = min(CT90);   [~,bestIdx_TW] = min(TW90);
[~,bestIdx_ER] = min(ER90);   specIdx = 1;

%% MIN–MAX RANGES ---------------------------------------------------------
fMin_CT = min(cdf_CT);  fMax_CT = max(cdf_CT);
fMin_TW = min(cdf_TW);  fMax_TW = max(cdf_TW);
fMin_ER = min(cdf_ER);  fMax_ER = max(cdf_ER);

%% PLOTTING ---------------------------------------------------------------
plotMetric(x_CT,cdf_CT,cdf_CT_dir,fMin_CT,fMax_CT,'$-\mathrm{Continuity}$', ...
           bestIdx_CT,specIdx,'b',[0 0 1],CTnegCell,'NegCT',figFolder,fileTag,nCodes);
plotMetric(x_TW,cdf_TW,cdf_TW_dir,fMin_TW,fMax_TW,'$-\mathrm{Trustworthiness}$', ...
           bestIdx_TW,specIdx,'r',[1 0 0],TWnegCell,'NegTW',figFolder,fileTag,nCodes);
plotMetric(x_ER,cdf_ER,cdf_ER_dir,fMin_ER,fMax_ER,'Positioning Error (m)', ...
           bestIdx_ER,specIdx,'g',[0 0.45 0],ERRcell,'PosErr',figFolder,fileTag,nCodes);

%% SAVE -------------------------------------------------------------------
bestIdxStruct = struct('bestCT',bestIdx_CT,'bestTW',bestIdx_TW,'bestER',bestIdx_ER);
save(fullfile(outFolder, ['MergedData_' fileTag '.mat']), ...
     'CTnegCell','TWnegCell','ERRcell', ...
     'cdf_CT','cdf_TW','cdf_ER', ...
     'bestIdxStruct','codebook_all','phi_options_1','phi_options_2','-v7.3');

fprintf('Finished. Data in  %s\n',outFolder);

%% ------------ LOCAL FUNCTION ------------------------------------------------
function plotMetric(xVals,cMat,cBase,fMin,fMax,labelStr,bestIdx,specIdx, ...
                    colorFill,colorBest,avgCell,filePrefix,figFolder,fileTag,nCodes)
    figure('Position',[100 100 800 600]); hold on; grid on; box on;
    hFill = fill([xVals fliplr(xVals)],[fMin fliplr(fMax)], ...
                 colorFill,'FaceAlpha',0.15,'EdgeColor','none');

    hDash = plot(nan,nan,'--','LineWidth',0.5,'Color',[.6 .6 .6 .6], ...
                 'DisplayName','All Other Codewords');
    for i=1:nCodes
        if any(i==[specIdx bestIdx]), continue; end
        plot(xVals,cMat(i,:),'--','LineWidth',0.5,'Color',[.6 .6 .6 .6],'HandleVisibility','off');
    end

    hSpec    = plot(xVals,cMat(specIdx,:),'k-.','LineWidth',2,  'DisplayName','Specular (\phi=0,0)');
    hBest    = plot(xVals,cMat(bestIdx,:),'-', 'LineWidth',2.5,'Color',colorBest,'DisplayName','Best Codeword');
    hBaseline= plot(xVals,cBase,          'k-','LineWidth',2.5,'DisplayName','No EMS (direct only)');

    specAvg = mean(avgCell{specIdx});  bestAvg = mean(avgCell{bestIdx});

    xlabel(labelStr,'Interpreter','latex','FontSize',18);
    ylabel('CDF','Interpreter','latex','FontSize',18);
    title(['Double RIS with AE: ' labelStr],'Interpreter','latex','FontSize',16);
    legend([hDash hFill hSpec hBest hBaseline], ...
        {'All Other Codewords',['Range of ' labelStr],'Specular','Best Codeword','No EMS (direct only)'}, ...
        'Interpreter','latex','FontSize',14,'Location','best');
    set(gca,'FontSize',14); ylim([0 1]);
    if contains(labelStr,'Error'), xlim([0 150]); else, xlim([min(xVals) 0]); end
    text(mean(xlim),0.05,sprintf('\\bf Spec avg=%.3f   Best avg=%.3f',specAvg,bestAvg), ...
        'Interpreter','latex','FontSize',14,'HorizontalAlignment','center');

    savefig([figFolder filePrefix '_' fileTag '.fig']);
    saveas(gcf,[figFolder filePrefix '_' fileTag '.png']);
    % % % close;
end

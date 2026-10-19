function params = init_parameters()
% INIT_PARAMETERS
% Centralized configuration for the DOUBLE EMS channel charting experiment.

    %% --------------------------------------------------------------------
    % REPRODUCIBILITY
    % ---------------------------------------------------------------------
    params.random_seed = 42;

    %% --------------------------------------------------------------------
    % PROJECT ROOT (NEW — IMPORTANT)
    % ---------------------------------------------------------------------
    % Automatically detect project root (parent of /src)
    thisFile = mfilename('fullpath');
    thisDir  = fileparts(thisFile);
    params.project_root = fileparts(thisDir);

    %% --------------------------------------------------------------------
    % DATASET
    % ---------------------------------------------------------------------
    params.dataPath = fullfile(params.project_root, ...
                               'data', 'HEMS_dataset.mat');

    %% --------------------------------------------------------------------
    % USER / EXPERIMENT CHOICES
    % ---------------------------------------------------------------------
    params.dissim_type = 'LE';
    params.Total_samples_used = 3200;
    params.Supervision_percentage = 0.15;
    params.supervised_Len = round( ...
        params.Total_samples_used * params.Supervision_percentage);

    params.WEIGHTING_METRIC = 'TW';
    params.SCENARIO_AGGREGATION = 'p90';

    %% --------------------------------------------------------------------
    % CHART TARGET
    % ---------------------------------------------------------------------
    % 'position' : original experiment, chart target = user (x,y) [m]
    % 'aoa'      : chart target = azimuth of the UE seen from EMS1 and EMS2
    %              [deg] (set by experiments/Main_AoA_Estimation.m)
    params.target = 'position';

    %% --------------------------------------------------------------------
    % EMS GEOMETRY FOR AoA LABELS (used only when params.target = 'aoa')
    % ---------------------------------------------------------------------
    % Verified from the dataset's channel phases:
    %   - both EMSs lie in walls perpendicular to x; the codebook phase
    %     gradient (element index m) runs along +y, the other axis along z;
    %   - EMS1 (x = -28.85) faces -x, EMS2 (x = -71.9) faces +x, i.e. both
    %     face the street between them, where the BS also is;
    %   - UE height 1.5 m (element spacing lambda/2).
    % EMS positions themselves are read from the dataset.
    params.UE_height  = 1.5;
    params.EMS_normal = [-1 0 0;  1 0 0];   % outward normal of EMS1, EMS2
    params.EMS_axis_m = [ 0 1 0;  0 1 0];   % codebook-gradient axis of EMS1, EMS2

    %% --------------------------------------------------------------------
    % EMS CONFIGURATION
    % ---------------------------------------------------------------------
    params.M = 60;
    params.N = 60;

    % % % % % % params.phi_options = deg2rad(0:4:40);
    params.phi_options = deg2rad(0:4:4);
    params.nPhasesPerEMS = numel(params.phi_options);

    %% --------------------------------------------------------------------
    % CHANNEL / SYSTEM PARAMETERS
    % ---------------------------------------------------------------------
    params.Par.Tx.ArrSize = 32;
    params.N_real = 500;

    %% --------------------------------------------------------------------
    % METRIC PARAMETERS
    % ---------------------------------------------------------------------
    params.K_fraction = 0.05;

    %% --------------------------------------------------------------------
    % EMBEDDING (t-SNE)
    % ---------------------------------------------------------------------
    params.tsne_perplexity = 150;

    %% --------------------------------------------------------------------
    % OUTPUT PATHS (UPDATED — now ROOT-BASED)
    % ---------------------------------------------------------------------
    params.output_base_fig = fullfile(params.project_root, 'Figures');
    params.output_base_res = fullfile(params.project_root, ...
                                      'results_double_EMS_cttw');

end
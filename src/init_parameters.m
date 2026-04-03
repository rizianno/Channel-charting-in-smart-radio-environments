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
    % EMS CONFIGURATION
    % ---------------------------------------------------------------------
    params.M = 60;
    params.N = 60;

    % % % params.phi_options = deg2rad(0:4:40);
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
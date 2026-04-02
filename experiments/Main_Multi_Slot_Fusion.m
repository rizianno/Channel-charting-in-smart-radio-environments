%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN SCRIPT — Modular version (1:1 with original behavior)
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

dataPath = fullfile(root, 'data', 'HEMS_dataset.mat');
figFolder = fullfile(root, 'Figures');
resultsFolder = fullfile(root, 'results');

Data  = load_dataset(params);
state = prepare_data(Data, params);

%% ------------------------------------------------------------------------
% CODEBOOK
%% ------------------------------------------------------------------------
codebook = generate_codebook(params);

fprintf('Number of codebooks (double EMS configurations): %d\n', codebook.nCodebooks);
fprintf('Weighting metric: %s (aggregation: %s)\n', ...
    params.WEIGHTING_METRIC, params.SCENARIO_AGGREGATION);

%% ------------------------------------------------------------------------
% PHASE 1 — CHANNEL CHARTING
%% ------------------------------------------------------------------------
results_phase1 = run_channel_charting(state, params, codebook);

fprintf('\nChannel charting completed for all %d codebooks.\n', ...
    codebook.nCodebooks);

%% ------------------------------------------------------------------------
% PHASE 2 — GLOBAL WEIGHTING (Scenario-level)
%% ------------------------------------------------------------------------
results_global = global_weighting(results_phase1, state, params);

%% ------------------------------------------------------------------------
% PHASE 2 — LOCAL WEIGHTING (Point-wise)
%% ------------------------------------------------------------------------
results_local = local_weighting(results_phase1, state, params);

%% ------------------------------------------------------------------------
% BEST CODEBOOK
%% ------------------------------------------------------------------------
best = evaluate_best_codebook(results_phase1, state, params);

% IMPORTANT: compute best error EXACTLY like original
err_best = sqrt(sum((best.xy - state.R_xy).^2, 2));
best.err = err_best(params.supervised_Len+1:end);

%% ------------------------------------------------------------------------
% METRICS (THIS IS YOUR FUNCTION)
%% ------------------------------------------------------------------------
metrics = compute_all_metrics(results_global, results_local, best, state, params);


%% ------------------------------------------------------------------------
% VISUALIZATION (IDENTICAL TO ORIGINAL)
%% ------------------------------------------------------------------------
visualize_results( ...
    results_phase1, ...
    results_global, ...
    results_local, ...
    best, ...
    metrics, ...
    state, ...
    params, ...
    codebook);

%% ------------------------------------------------------------------------
% SAVE RESULTS (IDENTICAL TO ORIGINAL)
%% ------------------------------------------------------------------------
save_results( ...
    results_phase1, ...
    results_global, ...
    results_local, ...
    best, ...
    metrics, ...
    state, ...
    params, ...
    codebook);
%% ------------------------------------------------------------------------
% PRINT SUMMARY (MATCH ORIGINAL)
%% ------------------------------------------------------------------------
fprintf('\n=== SUMMARY TABLE ===\n');
fprintf('---------------------------------------------------------------------------------\n');
fprintf('Metric                  Best Codebook  Global Weighted  Local Weighted\n');
fprintf('---------------------------------------------------------------------------------\n');

fprintf('Mean Error (m)          %8.2f       %8.2f           %8.2f\n', ...
    best.error, ...
    mean(results_global.err), ...
    mean(results_local.err));

fprintf('Median Error (m)        %8.2f       %8.2f           %8.2f\n', ...
    median(best.err), ...
    median(results_global.err), ...
    median(results_local.err));

fprintf('90th Percentile (m)     %8.2f       %8.2f           %8.2f\n', ...
    prctile(best.err,90), ...
    prctile(results_global.err,90), ...
    prctile(results_local.err,90));

fprintf('Mean CT                 %8.3f       %8.3f           %8.3f\n', ...
    mean(metrics.CT_best), ...
    mean(metrics.CT_global), ...
    mean(metrics.CT_local));

fprintf('Mean TW                 %8.3f       %8.3f           %8.3f\n', ...
    mean(metrics.TW_best), ...
    mean(metrics.TW_global), ...
    mean(metrics.TW_local));

fprintf('---------------------------------------------------------------------------------\n');

%% ------------------------------------------------------------------------
% IMPROVEMENT (same as original logic)
%% ------------------------------------------------------------------------
improve_global = (best.error - mean(results_global.err)) / best.error * 100;
improve_local  = (best.error - mean(results_local.err)) / best.error * 100;

fprintf('\nImprovement vs Best Single Codebook:\n');
fprintf('  Global weighted:  %.1f%%\n', improve_global);
fprintf('  Local weighted:   %.1f%%\n', improve_local);

fprintf('\n=== DONE ===\n');





function visualize_results(results_phase1, results_global, results_local, best, metrics, state, params, codebook)
% VISUALIZE_RESULTS
% Fully identical visualization pipeline

    fprintf('\n=== Creating Visualizations ===\n');

    Len = state.Len;
    nCodebooks = size(results_phase1.PosEstimates,3);

    %% Reconstruct commonly used variables
    ErrorPerCodebook = results_phase1.ErrorPerCodebook;
    ErrorPerUser     = results_phase1.ErrorPerUser;

    R_xy = state.R_xy;

    err_best = best.err;
    err_global = results_global.err;
    err_local  = results_local.err;

    %% --------------------------------------------------------------------
    % FIGURE 1: METRICS PER CODEBOOK
    % --------------------------------------------------------------------
    figure('Position', [100, 100, 1400, 900]);

    subplot(2,2,1);
    plot(1:nCodebooks, ErrorPerCodebook, '-', 'Color', [0.7 0.7 0.7]); hold on;
    plot(best.idx, best.error, 'bo', 'LineWidth', 2);

    yline(mean(err_global), 'r--', 'LineWidth', 2);
    yline(mean(err_local),  'g--', 'LineWidth', 2);
    yline(best.error,       'b--', 'LineWidth', 2);

    xlabel('Codebook Index');
    ylabel('Average Error (m)');
    title('Positioning Error per Codebook');
    grid on;

    subplot(2,2,2);
    scatter3(codebook.all(:,1)*180/pi, codebook.all(:,2)*180/pi, ...
        ErrorPerCodebook, 50, ErrorPerCodebook, 'filled');
    colorbar; colormap('hot'); grid on;

    subplot(2,2,3);
    yyaxis left
    plot(results_phase1.CT_PerCodebook_mean);
    yyaxis right
    plot(results_phase1.TW_PerCodebook_mean);
    title('CT/TW per Codebook'); grid on;

    subplot(2,2,4);
    histogram(results_global.weights,30);
    title('Global Weight Distribution');

    %% --------------------------------------------------------------------
    % FIGURE 2: CDF
    % --------------------------------------------------------------------
    figure; hold on; grid on;

    cb_subsample = 1:max(1,floor(nCodebooks/50)):nCodebooks;

    for cb = cb_subsample
        err_cb = ErrorPerUser(params.supervised_Len+1:end, cb);
        plot(sort(err_cb), linspace(0,1,length(err_cb)), ...
            'Color',[0.85 0.85 0.85]);
    end

    plot(sort(err_best), linspace(0,1,length(err_best)), 'b','LineWidth',2);
    plot(sort(err_global), linspace(0,1,length(err_global)), 'r','LineWidth',2);
    plot(sort(err_local), linspace(0,1,length(err_local)), 'g','LineWidth',2);

    yline(0.9,'k--');

    xlabel('Error (m)');
    ylabel('CDF');
    title('Error Distribution');
end



function save_results(results_phase1, results_global, results_local, best, metrics, state, params, codebook)
% SAVE_RESULTS
% Stores everything in same structure as original script

    fprintf('\n=== Saving Results ===\n');

    results = struct();

    %% Phase 1
    results.codebook_all = codebook.all;
    results.phi_options  = params.phi_options;

    results.ErrorPerCodebook = results_phase1.ErrorPerCodebook;
    results.ErrorPerUser     = results_phase1.ErrorPerUser;
    results.PosEstimates     = results_phase1.PosEstimates;

    %% Global
    results.global_weights = results_global.weights;
    results.xy_global = results_global.xy;
    results.err_global = results_global.err;
    results.CT_global = metrics.CT_global;
    results.TW_global = metrics.TW_global;

    %% Local
    results.local_weights = results_local.weights;
    results.xy_local = results_local.xy;
    results.err_local = results_local.err;
    results.CT_local = metrics.CT_local;
    results.TW_local = metrics.TW_local;

    %% Best
    results.best_codebook_idx = best.idx;
    results.xy_best = best.xy;
    results.err_best = best.err;
    results.CT_best = metrics.CT_best;
    results.TW_best = metrics.TW_best;

    %% Meta
    results.R_xy = state.R_xy;
    results.supervised_Len = params.supervised_Len;
    results.weighting_metric = params.WEIGHTING_METRIC;
    results.scenario_aggregation = params.SCENARIO_AGGREGATION;

    %% Save
    if ~exist('results','dir'), mkdir('results'); end
    save('results/results.mat','results','-v7.3');

    fprintf('Results saved.\n');

end
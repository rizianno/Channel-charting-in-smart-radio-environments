function results_global = global_weighting(results_phase1, state, params)
% GLOBAL_WEIGHTING
% Performs global (scenario-level) weighted fusion of embeddings.
%
% Each EMS configuration is assigned a single global weight based on
% aggregated CT/TW metrics across all unlabeled users.
%
% This corresponds to:
%   Sec. IV-D2 (global fusion strategy)

    fprintf('\n=== PHASE 2A: Global Weighted Averaging ===\n');

    Len = state.Len;
    nCodebooks = size(results_phase1.PosEstimates,3);

    %% --------------------------------------------------------------------
    % SELECT GLOBAL QUALITY METRIC
    % ---------------------------------------------------------------------
    if strcmpi(params.WEIGHTING_METRIC, 'TW')
        if strcmpi(params.SCENARIO_AGGREGATION, 'mean')
            global_quality = results_phase1.TW_PerCodebook_mean;
        else
            global_quality = results_phase1.TW_PerCodebook_p90;
        end
    else  % CT
        if strcmpi(params.SCENARIO_AGGREGATION, 'mean')
            global_quality = results_phase1.CT_PerCodebook_mean;
        else
            global_quality = results_phase1.CT_PerCodebook_p90;
        end
    end

    %% --------------------------------------------------------------------
    % NORMALIZE WEIGHTS
    % ---------------------------------------------------------------------
    global_weights = global_quality / sum(global_quality);

    %% --------------------------------------------------------------------
    % LOGGING
    % ---------------------------------------------------------------------
    fprintf('Global weights (based on %s %s):\n', ...
        params.WEIGHTING_METRIC, params.SCENARIO_AGGREGATION);

    fprintf('  Min weight: %.6f, Max weight: %.6f\n', ...
        min(global_weights), max(global_weights));

    [~, sorted_idx] = sort(global_weights, 'descend');
    fprintf('  Top 5 codebooks:\n');
    for i = 1:min(5,nCodebooks)
        cb = sorted_idx(i);
        fprintf('    Codebook %3d: weight = %.6f (quality = %.4f)\n', ...
            cb, global_weights(cb), global_quality(cb));
    end

    %% --------------------------------------------------------------------
    % WEIGHTED FUSION
    % ---------------------------------------------------------------------
    xy_weighted = zeros(Len,2);

    for uu = 1:Len
        for cb = 1:nCodebooks
            xy_weighted(uu,:) = xy_weighted(uu,:) + ...
                global_weights(cb) * results_phase1.PosEstimates(uu,:,cb);
        end
    end

    % Clamp labeled users
    xy_weighted(1:params.supervised_Len,:) = ...
        state.R_xy(1:params.supervised_Len,:);

    %% --------------------------------------------------------------------
    % ERROR EVALUATION
    % ---------------------------------------------------------------------
    err = sqrt(sum((xy_weighted - state.R_xy).^2,2));
    err_unlabeled = err(params.supervised_Len+1:end);

    fprintf('\nGlobal Weighted Results:\n');
    fprintf('  Mean error:        %.2f m\n', mean(err_unlabeled));
    fprintf('  Median error:      %.2f m\n', median(err_unlabeled));
    fprintf('  90th percentile:   %.2f m\n', prctile(err_unlabeled,90));

    %% --------------------------------------------------------------------
    % STORE OUTPUT
    % ---------------------------------------------------------------------
    results_global.weights = global_weights;
    results_global.xy = xy_weighted;
    results_global.err = err_unlabeled;

end
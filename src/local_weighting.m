function results_local = local_weighting(results_phase1, state, params)
% LOCAL_WEIGHTING
% Performs local (per-user) weighted fusion of embeddings.
%
% Each user assigns weights to EMS configurations based on its own
% CT/TW metric values.
%
% This corresponds to:
%   Sec. IV-D2 (local fusion strategy)

    fprintf('\n=== PHASE 2B: Local Weighted Averaging ===\n');

    Len = state.Len;
    nCodebooks = size(results_phase1.PosEstimates,3);

    %% --------------------------------------------------------------------
    % SELECT LOCAL QUALITY METRIC
    % ---------------------------------------------------------------------
    if strcmpi(params.WEIGHTING_METRIC, 'TW')
        local_quality = results_phase1.TW_PerUser;
    else
        local_quality = results_phase1.CT_PerUser;
    end

    %% --------------------------------------------------------------------
    % COMPUTE PER-USER WEIGHTS
    % ---------------------------------------------------------------------
    local_weights = zeros(Len, nCodebooks);

    for uu = 1:Len

        if uu > params.supervised_Len  % unlabeled users

            quality_vec = local_quality(uu,:);

            % Avoid zero division / numerical issues
            quality_vec = max(quality_vec, 1e-6);

            local_weights(uu,:) = quality_vec / sum(quality_vec);

        else
            % Labeled users (will be clamped anyway)
            local_weights(uu,:) = 1/nCodebooks;
        end
    end

    %% --------------------------------------------------------------------
    % WEIGHTED FUSION
    % ---------------------------------------------------------------------
    xy_weighted = zeros(Len,2);

    for uu = 1:Len
        for cb = 1:nCodebooks
            xy_weighted(uu,:) = xy_weighted(uu,:) + ...
                local_weights(uu,cb) * results_phase1.PosEstimates(uu,:,cb);
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

    fprintf('\nLocal Weighted Results:\n');
    fprintf('  Mean error:        %.2f m\n', mean(err_unlabeled));
    fprintf('  Median error:      %.2f m\n', median(err_unlabeled));
    fprintf('  90th percentile:   %.2f m\n', prctile(err_unlabeled,90));

    %% --------------------------------------------------------------------
    % STORE OUTPUT
    % ---------------------------------------------------------------------
    results_local.weights = local_weights;
    results_local.xy = xy_weighted;
    results_local.err = err_unlabeled;

end
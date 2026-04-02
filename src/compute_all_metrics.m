function metrics = compute_all_metrics(results_global, results_local, best, state, params)
% COMPUTE_ALL_METRICS
% Computes CT/TW metrics for:
%   - global weighting
%   - local weighting
%   - best single EMS configuration

    fprintf('\n=== Computing CT/TW Metrics ===\n');

    Len = state.Len;

    KK = 0.05 * Len;
    DistT = squareform(pdist(state.R_xy));

    %% GLOBAL
    DistZ_global = squareform(pdist(results_global.xy));
    TC_global = drscore_pointwise(DistT, DistZ_global, KK, 'TC');

    metrics.CT_global = TC_global(params.supervised_Len+1:end,1);
    metrics.TW_global = TC_global(params.supervised_Len+1:end,2);

    %% LOCAL
    DistZ_local = squareform(pdist(results_local.xy));
    TC_local = drscore_pointwise(DistT, DistZ_local, KK, 'TC');

    metrics.CT_local = TC_local(params.supervised_Len+1:end,1);
    metrics.TW_local = TC_local(params.supervised_Len+1:end,2);

    %% BEST CODEBOOK
    DistZ_best = squareform(pdist(best.xy));
    TC_best = drscore_pointwise(DistT, DistZ_best, KK, 'TC');

    metrics.CT_best = TC_best(params.supervised_Len+1:end,1);
    metrics.TW_best = TC_best(params.supervised_Len+1:end,2);

end
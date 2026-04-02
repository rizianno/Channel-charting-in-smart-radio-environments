function best = evaluate_best_codebook(results_phase1, state, params)
% EVALUATE_BEST_CODEBOOK
% Identifies the single best EMS configuration based on mean error.

    [best.error, best.idx] = min(results_phase1.ErrorPerCodebook);

    best.xy = squeeze(results_phase1.PosEstimates(:,:,best.idx));

    % Compute error explicitly (for consistency with others)
    err = sqrt(sum((best.xy - state.R_xy).^2,2));
    best.err = err(params.supervised_Len+1:end);

    fprintf('\nBest Single Codebook:\n');
    fprintf('  Index: %d\n', best.idx);
    fprintf('  Mean error: %.2f m\n', best.error);
    fprintf('  90th percentile: %.2f m\n', prctile(best.err,90));

end
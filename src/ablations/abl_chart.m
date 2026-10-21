function out = abl_chart(D, labels, nA, valid, perplexity, Kfrac)
% ABL_CHART  One semi-supervised t-SNE chart (tsne_sem, unchanged) plus the per-user quality
% used by the fusion.
%   D      : [n x n] input dissimilarities, users ordered with the nA anchors first
%   labels : [n x d] target (azimuths); only the anchors' rows are used by t-SNE
%   valid  : [n x 1] users whose label is defined (for the ground-truth quality)
% Quality, as in run_channel_charting + global/local_weighting (WEIGHTING_METRIC = 'TW'):
% column 2 of drscore_pointwise, which is the continuity (the original code stores it under
% the name TW).
%   q_orig : vs the TRUE target distances of all users (uses the unlabeled users' labels,
%            as the original code does); NaN where the label is undefined
%   q_free : vs the input dissimilarities D (label-free alternative)
    n   = size(D,1);
    est = tsne_sem(D, labels, 1:nA, perplexity);
    Dz  = squareform(pdist(est));

    vi = find(valid);
    q_orig = nan(n,1);
    TC = abl_drscore(squareform(pdist(labels(vi,:))), Dz(vi,vi), Kfrac*numel(vi));
    q_orig(vi) = TC(:,2);

    TC2 = abl_drscore(D, Dz, Kfrac*n);
    out.est    = est;
    out.q_orig = q_orig;
    out.q_free = TC2(:,2);
end

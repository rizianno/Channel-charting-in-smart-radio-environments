function out = srv_chart_twct(D, labels, nA, perplexity, Kfrac)
% SRV_CHART_TWCT  One semi-supervised t-SNE chart (tsne_sem, unchanged) and its label-free
% per-user trustworthiness and continuity: neighbourhoods in the input dissimilarities D vs in
% the chart (drscore_pointwise convention: column 1 = trustworthiness, column 2 = continuity).
%   D      : [n x n] input dissimilarities, the nA anchors first
%   labels : [n x d] targets (only the anchors' rows are used by t-SNE)
    n = size(D,1);
    out.est = tsne_sem(D, labels, 1:nA, perplexity);
    TC = abl_drscore(D, squareform(pdist(out.est)), Kfrac*n);
    out.tw = TC(:,1);
    out.ct = TC(:,2);
end

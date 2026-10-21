function [g, l] = abl_fuse(EST, Q, nA, lab)
% ABL_FUSE  Global and local quality-weighted fusion of per-sweep charts, as global_weighting
% (p90 of the unlabeled users' quality per chart) and local_weighting (per-user quality, floored
% at 1e-6). Users without a quality value get the global weights. Anchors keep their labels.
%   EST : [n x d x C] chart estimates,  Q : [n x C] per-user quality,  lab : [n x d]
    [n, ~, C] = size(EST);
    unl = nA+1:n;
    gq = zeros(1,C);
    for c = 1:C
        qc = Q(unl,c);
        gq(c) = prctile(qc(~isnan(qc)), 90);
    end
    wg = gq / sum(gq);
    g  = sum(EST .* reshape(wg,1,1,C), 3);

    Wl = max(Q, 1e-6);
    Wl = Wl ./ sum(Wl, 2);
    bad = any(isnan(Wl), 2);
    Wl(bad,:) = repmat(wg, nnz(bad), 1);
    l  = sum(EST .* reshape(Wl,n,1,C), 3);

    g(1:nA,:) = lab(1:nA,:);
    l(1:nA,:) = lab(1:nA,:);
end

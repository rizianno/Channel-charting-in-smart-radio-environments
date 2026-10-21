function C = abl_charts(U, structure)
% ABL_CHARTS  t-SNE charts that make up one user/chart structure.
%   .id, .users (dataset indices, anchors first), .nA, .cols (azimuth columns = labels),
%   .labels, .valid (label defined), .perEMS (true for 1-D per-surface charts)
    switch structure
        case 'S1'      % joint (az1, az2), users in front of both
            C = mk('S1', U.anchorsJoint, find(U.both), [1 2], U.both, U);
        case 'S2'      % one 1-D chart per surface, users in front of that surface
            C = [mk('S2e1', U.anchorsEMS{1}, find(U.front(:,1)), 1, U.front(:,1), U), ...
                 mk('S2e2', U.anchorsEMS{2}, find(U.front(:,2)), 2, U.front(:,2), U)];
        case 'S3'      % joint, all users in t-SNE, same anchors as S1
            C = mk('S3', U.anchorsJoint, (1:U.n).', [1 2], U.both, U);
        case 'S4'      % one 1-D chart per surface on the S1 users and anchors
            C = [mk('S4e1', U.anchorsJoint, find(U.both), 1, U.both, U), ...
                 mk('S4e2', U.anchorsJoint, find(U.both), 2, U.both, U)];
        otherwise
            error('Unknown structure %s', structure);
    end
end

function c = mk(id, anchors, members, cols, validAll, U)
    c.id     = id;
    c.users  = [anchors(:); setdiff(members(:), anchors(:), 'stable')];
    c.nA     = numel(anchors);
    c.cols   = cols;
    c.valid  = validAll(c.users);
    c.labels = U.az(c.users, cols);
    c.labels(~c.valid,:) = NaN;
    c.perEMS = isscalar(cols);
end

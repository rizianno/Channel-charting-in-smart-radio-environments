function J = abl_jobs(scn, K, charts, cat, noiseName)
% ABL_JOBS  t-SNE jobs needed by one scenario at codebook size K x K.
% F0: one job per (chart, configuration), shared by every K that contains the configuration.
% F1-F3: one job per (chart, K).
%   .key, .chart, .feature, .direct, .phi [S x 2] configurations (deg),
%   .off [S x 2] OFF acquisition for each configuration (F3), .group [S x 1] centring groups (F2)
    lev = abl_levels(cat, scn, K);
    [p1, p2] = ndgrid(lev, lev);
    P = [p1(:) p2(:)];                    % S = K^2 configurations
    J = struct('key',{}, 'chart',{}, 'feature',{}, 'direct',{}, 'phi',{}, 'off',{}, 'group',{});
    for c = charts
        if strcmp(scn.feature, 'F0')
            for s = 1:size(P,1)
                J(end+1) = mk(sprintf('%s|%s|F0|%s', noiseName, c.id, ...
                              abl_stkey(scn.direct, P(s,:))), c, 'F0', scn.direct, P(s,:)); %#ok<AGROW>
            end
        else
            cbk = '';  if isfield(scn,'codebook') && ~strcmp(scn.codebook,'narrow'), cbk = ['|' scn.codebook]; end
            j = mk(sprintf('%s|%s|%s|d%d|K%d%s', noiseName, c.id, scn.feature, scn.direct, K, cbk), ...
                   c, scn.feature, scn.direct, P);
            if c.perEMS
                e = c.cols;  o = 3 - e;
                j.off = P;  j.off(:,e) = NaN;               % switch off only this surface
                [~,~,j.group] = unique(P(:,o));            % centre along this surface's sweep
            else
                j.off = NaN(size(P));                      % switch off both surfaces
                j.group = ones(size(P,1),1);               % centre over all configurations
            end
            J(end+1) = j; %#ok<AGROW>
        end
    end
end

function j = mk(key, c, f, d, P)
    j = struct('key',key, 'chart',c, 'feature',f, 'direct',d, 'phi',P, 'off',[], 'group',[]);
end

function [trajIdx, trajXY, speed, accel, timeVec] = ...
         buildTrajectory_smooth(R_xy, totalTime, dt, v_max, a_max, startIdx, K)
% BUILDTRAJECTORY_SMOOTH
% Generates a long, smooth random walk on a 2-D measurement grid.
%
% The walker moves between grid points (nodes) subject to:
%   • a speed cap        v_max  [m/s]  → maximum single-step distance d_max
%   • an acceleration cap a_max [m/s²] → Δv per tick is bounded
%   • a heading-alignment score        → neighbours that continue the current
%     heading direction are preferred (cosine-similarity greedy selection)
%
% At each tick the walker picks the reachable neighbour with the highest
% cosine similarity to the previous heading while satisfying the
% acceleration constraint.  If no neighbour satisfies the acceleration
% bound, the constraint is relaxed and the nearest neighbour is chosen.
%
% Inputs:
%   R_xy       : [N x 2]   grid node coordinates [m]
%   totalTime  : scalar    total simulation time [s]
%   dt         : scalar    time step [s]
%   v_max      : scalar    maximum speed [m/s]
%                          If empty, auto-selected as 1.2 × median(nearest-
%                          neighbour distance) / dt
%   a_max      : scalar    acceleration bound [m/s²]
%   startIdx   : integer   starting node index (1 … N)
%   K          : integer   (optional) k-NN fallback for connectivity;
%                          default = 4
%
% Outputs:
%   trajIdx  : [nTicks x 1]  node indices visited at each tick
%   trajXY   : [nTicks x 2]  corresponding (x,y) positions [m]
%   speed    : [nTicks x 1]  instantaneous speed at each tick [m/s]
%   accel    : [nTicks x 1]  instantaneous acceleration [m/s²]
%   timeVec  : [nTicks x 1]  time vector 0 : dt : totalTime  [s]

    if nargin < 7, K = 4; end

    %% ------------------------------------------------------------------
    % Auto-select v_max if not provided
    % ------------------------------------------------------------------
    if isempty(v_max)
        nearDist = min(pdist2(R_xy, R_xy) + diag(inf(size(R_xy,1),1)), [], 2);
        v_max    = 1.2 * median(nearDist) / dt;
        fprintf('buildTrajectory_smooth: auto v_max = %.3f m/s  (d_max = %.2f m)\n', ...
                v_max, v_max*dt);
    end

    nTicks   = floor(totalTime / dt) + 1;
    d_max    = v_max * dt;
    a_max_dt = a_max * dt;             % max speed change per tick [m/s]

    %% ------------------------------------------------------------------
    % Build radius-graph adjacency, with K-NN fallback for isolated nodes
    % ------------------------------------------------------------------
    D   = squareform(pdist(R_xy));
    adj = (D > 0) & (D <= d_max);

    if any(sum(adj, 2) == 0)
        [~, ord] = sort(D, 2);
        for i = 1:size(D, 1)
            adj(i, ord(i, 2:K+1)) = true;
        end
    end

    %% ------------------------------------------------------------------
    % Allocate outputs
    % ------------------------------------------------------------------
    trajIdx = zeros(nTicks, 1);
    trajIdx(1) = startIdx;
    speed   = zeros(nTicks, 1);
    accel   = zeros(nTicks, 1);
    prevDir = [];                      % heading unit-vector from last step

    %% ------------------------------------------------------------------
    % Main walk
    % ------------------------------------------------------------------
    for t = 2:nTicks

        cur  = trajIdx(t-1);
        nbrs = find(adj(cur, :));

        if isempty(nbrs)
            [~, nbrs] = min(D(cur, :));
        end

        v_prev    = (t > 2) * speed(t-1);
        bestScore = -inf;
        chosen    = [];

        % Greedy: pick the neighbour that (a) respects accel bound and
        % (b) maximises cosine similarity with the current heading
        for n = nbrs(randperm(numel(nbrs)))
            dist = D(cur, n);
            v    = dist / dt;
            if abs(v - v_prev) > a_max_dt + 1e-12, continue; end

            if isempty(prevDir)
                score = 1;             % first step: no heading preference
            else
                dirVec = (R_xy(n,:) - R_xy(cur,:)) / dist;
                score  = dirVec * prevDir';
            end

            if score > bestScore
                bestScore = score;
                chosen    = [n, dist, v];
            end
        end

        if isempty(chosen)             % relax acceleration constraint
            n    = nbrs(1);
            dist = D(cur, n);
            v    = dist / dt;
        else
            n    = chosen(1);
            dist = chosen(2);
            v    = chosen(3);
        end

        trajIdx(t) = n;
        speed(t)   = v;
        accel(t)   = (v - v_prev) / dt;
        prevDir    = (R_xy(n,:) - R_xy(cur,:)) / dist;
    end

    trajXY  = R_xy(trajIdx, :);
    timeVec = (0:nTicks-1).' * dt;

end

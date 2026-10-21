function TC = abl_drscore(D0, D, K)
% ABL_DRSCORE  Same as drscore_pointwise(D0, D, K, 'TC'), i.e. [Trustworthiness Continuity]
% per point, but with O(N) rank look-ups instead of find() (identical results, faster).
    N  = size(D0,1);
    k  = floor(K);
    TW = zeros(N,1); CT = zeros(N,1);
    rV = zeros(1,N); rU = zeros(1,N);
    for i = 1:N
        d0 = D0(i,:); d0(i) = 1e6;
        d  = D(i,:);  d(i)  = 1e6;
        [~,Vn] = sort(d0); [~,Un] = sort(d);
        rV(Vn) = 1:N; rU(Un) = 1:N;
        Vk = Vn(1:k); Uk = Un(1:k);
        inU = false(1,N); inU(Uk) = true;
        inV = false(1,N); inV(Vk) = true;
        CT(i) = sum(rU(Vk(~inU(Vk))) - K);   % rank in the map of original neighbours that left
        TW(i) = sum(rV(Uk(~inV(Uk))) - K);   % rank in the original space of intruders in the map
    end
    c  = 2/(K*(2*N-3*K-1));
    TC = [1 - c*TW, 1 - c*CT];
end

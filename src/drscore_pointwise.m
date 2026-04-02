function [score,Con ]= drscore_pointwise(D0,D,K,option)
N = size(D0,1);
switch option
    case 'TC'
        Continuity = zeros(N,1);
        Trustworthiness = zeros(N,1);
        Con = zeros(N,1);
        for i = 1:N
        
            dist_v0 = D0(i,:);  dist_v0(i) = 1e6;
            dist_v  = D(i,:);   dist_v(i)  = 1e6;
            [~,Vn] = sort(dist_v0); Vk = Vn(1:K); % the set of points that are among the K nearest neighbors in the orignal space
            [~,Un] = sort(dist_v);  Uk = Un(1:K); % the set of points that are among the K nearest neighbors in the mapping space
            VVk = setdiff(Vk,Uk);
            UUk = setdiff(Uk,Vk);
            for index=1:length(VVk)
                j = VVk(index);
                rank_ij = find(Un==j);  % the rank of j in the mapping space to point i
                Continuity(i) = Continuity(i) + (rank_ij - K);
                Con(i) = Con(i) + (rank_ij - K);
            end
            for index=1:length(UUk)
                j = UUk(index);
                rank_ij = find(Vn==j);  % the rank of j in the orignal space to point i
                Trustworthiness(i) = Trustworthiness(i) + (rank_ij - K);
            end
        Continuity(i)      = 1-(2/(K*(2*N-3*K-1)))*Continuity(i);
        Trustworthiness(i) = 1-(2/(K*(2*N-3*K-1)))*Trustworthiness(i);
        end
        
        % Continuity      = 1-(2/(N*K*(2*N-3*K-1)))*Continuity;
        % Trustworthiness = 1-(2/(N*K*(2*N-3*K-1)))*Trustworthiness;
        
        score = [Trustworthiness Continuity];

    case 'KruskalStress'
        KS = zeros(N,1);
        Con = zeros(N,1);
        for i =1:N
            lamda  = sum(D0(i,:).*D(i,:))./sum(D(i,:).^2);
            KS(i)  = sqrt(sum((D0(i,:)-lamda.*D(i,:)).^2 )/sum(D0(i,:).^2));
        end
        score = KS;

    otherwise
        disp('other value')
end
end

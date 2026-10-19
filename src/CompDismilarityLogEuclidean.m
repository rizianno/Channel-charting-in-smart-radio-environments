function [D_LogEuc, Feat_sav]= CompDismilarityLogEuclidean(R_Cov_BS)
% compute the dissmialrity based on Log-Euclidean  distance
% comput dissimilarity between two features (here two users)
%
% D_LogEuc(i,j) = || logm(R_i + tol*I) - logm(R_j + tol*I) ||_F
%
% Same result as the original double loop, but logm is computed once per
% user (instead of once per pair) and the pairwise Frobenius distances are
% taken between the vectorized logs: for complex matrices,
% ||A - B||_F = || [real(vec A), imag(vec A)] - [real(vec B), imag(vec B)] ||_2.

Nu_UE    = size(R_Cov_BS,1);
N_Ant    = size(R_Cov_BS,2);

%  find indices of  the diagonal, lower diagonal (for the feature output)
n = 1:N_Ant;
n_d= N_Ant*(n-1)+n; % diagonal element
n_l=[];
for i=1:N_Ant
    n_i=N_Ant*(i-1)+i+1:N_Ant*(i-1)+N_Ant;
    n_l=[n_l, n_i];
end

esp_tol=1e-6;

L        = complex(zeros(Nu_UE, N_Ant*N_Ant));   % vec(logm(R_u)) per user
Feat_sav = zeros(Nu_UE, numel(n_d) + 2*numel(n_l));
for i_u = 1: Nu_UE
    R_iu = squeeze(R_Cov_BS(i_u,:,:));
    R_iu = R_iu + esp_tol*eye(N_Ant);
    R_iu_tangSpac = logm(R_iu);
    L(i_u,:) = R_iu_tangSpac(:).';
    Feat_sav(i_u,:) = [real(R_iu_tangSpac(n_d)), sqrt(2)*real(R_iu_tangSpac(n_l)), ...
                       sqrt(2)*imag(R_iu_tangSpac(n_l))];
end

D_LogEuc = squareform(pdist([real(L), imag(L)]));

end

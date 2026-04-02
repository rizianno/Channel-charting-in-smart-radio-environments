function [D_LogEuc, Feat_sav]= CompDismilarityLogEuclidean(R_Cov_BS)
% compute the dissmialrity based on Log-Euclidean  distance
% comput dissimilarity between two features (here two users)


Nu_UE    = size(R_Cov_BS,1);
N_Ant    = size(R_Cov_BS,2);
%beata= 4; % 4, 16
%  find indices of  the diagonal, lower diagonal
n = 1:N_Ant;
n_d= N_Ant*(n-1)+n; % diagonal element
n_l=[];
for i=1:N_Ant
    n_i=N_Ant*(i-1)+i+1:N_Ant*(i-1)+N_Ant;
    n_l=[n_l, n_i];
end

esp_tol=1e-6;
D_LogEuc  = zeros(Nu_UE,Nu_UE);
% D_CMD  = zeros(Nu_UE,Nu_UE);
Feat_sav=[];
parfor i_u = 1: Nu_UE
    disp(['calculating dissimilarity for user ',num2str(i_u), ' out of ',num2str(Nu_UE)])
    R_iu = squeeze(R_Cov_BS(i_u,:,:)); % size:  nBSbeams * nBSbeams
    %R_iu= R_iu/norm( R_iu,'fro'); % normalization by fro norm
    %R_iu = R_iu/trace(R_iu);       % normalization by trace of Cov matrix
    %R_iu= R_iu/norm( R_iu,'fro')^(1+1/(2*beata)); %with scaling
    R_iu=  R_iu+esp_tol*eye(N_Ant); % we do not need dissimilarity of users with themselves
    R_iu_tangSpac=logm(R_iu);
    %[eigVec_iu, eigvalu_iu]=eig(R_iu);
    %R_iu_tangSpac= eigVec_iu*diag(log(diag(eigvalu_iu)))*eigVec_iu';
    Feat_iu= [real(R_iu_tangSpac(n_d)), sqrt(2)*real(R_iu_tangSpac(n_l))...
        ,sqrt(2)*imag(R_iu_tangSpac(n_l))]';
    Feat_sav = [Feat_sav;Feat_iu'];
    for j_u = 1: Nu_UE
        R_ju = squeeze(R_Cov_BS(j_u,:,:));
        %R_ju= R_ju/norm( R_ju,'fro'); %withot normalization
        %R_ju= R_ju/norm( R_ju,'fro')^(1+1/(2*beata));  %with scaling
        R_ju=  R_ju+esp_tol*eye(N_Ant);
        %[eigVec_ju, eigvalu_ju]=eig(R_ju);
        %R_ju_tangSpac= eigVec_ju*diag(log(diag(eigvalu_ju)))*eigVec_ju';
        R_ju_tangSpac=logm(R_ju);
        D_LogEuc( i_u, j_u)   = norm(R_iu_tangSpac- R_ju_tangSpac,'fro');
        % D_LogEuc(j_u, i_u)  =  D_LogEuc(i_u, j_u);
        % D_CMD( i_u, j_u) = 1/2*norm(R_iu_tangSpac- R_ju_tangSpac,'fro');
        % D_CMD(j_u, i_u)  =  D_LogEuc(i_u, j_u);
        
        
    end
    if ~rem(i_u,100)
        disp(['considered users', num2str(i_u)])
    end
end


end


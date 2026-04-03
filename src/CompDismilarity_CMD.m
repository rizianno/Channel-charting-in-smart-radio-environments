function [D_CMD, Feat_sav] = CompDismilarity_CMD(R_Cov_BS)
% COMPDISMILARITY_CMD
% Computes pairwise dissimilarity using the Correlation Matrix Distance (CMD).
%
%   CMD(A, B) = 1 - trace(A * B) / ( ||A||_F * ||B||_F )
%
% For Hermitian positive semi-definite covariance matrices A and B,
% trace(A*B) is real and non-negative, so CMD lies in [0, 1].
%
% This function is a drop-in peer of CompDismilarityLogEuclidean:
%   same input format, same output format.
%
% Inputs:
%   R_Cov_BS  : [Nu_UE x N_Ant x N_Ant] array of covariance matrices
%               (complex Hermitian PSD, one per user)
%
% Outputs:
%   D_CMD     : [Nu_UE x Nu_UE] symmetric dissimilarity matrix,
%               values in [0, 1]
%   Feat_sav  : [] — placeholder for API compatibility with
%               CompDismilarityLogEuclidean (no tangent-space features
%               are extracted by CMD)
%
% Notes on the vectorised computation:
%   For Hermitian PSD A and B, let a = A(:) and b = B(:) be column-major
%   flat vectors.  Then:
%       trace(A * B) = sum_{r,c} A(r,c) * conj(B(r,c))
%                   = real( a.' * conj(b) )
%                   = real( Gram matrix entry (i,j) of  Mats * Mats' )
%   where Mats stores each flattened matrix as a row.
%   The Frobenius norm ||A||_F = sqrt( real( a.' * conj(a) ) )
%                              = sqrt( diagonal of Gram matrix ).

    Nu_UE = size(R_Cov_BS, 1);
    N_Ant = size(R_Cov_BS, 2);

    Feat_sav = [];

    %% ------------------------------------------------------------------
    % Flatten covariance matrices to rows of a [Nu_UE x N_Ant^2] matrix
    % ------------------------------------------------------------------
    Mats = reshape(R_Cov_BS, Nu_UE, N_Ant * N_Ant);   % complex

    %% ------------------------------------------------------------------
    % Gram matrix: GM(i,j) = trace(R_i * R_j)  (real for Hermitian PSD)
    % Diagonal:   GM(i,i) = ||R_i||_F^2
    % ------------------------------------------------------------------
    GM = real(Mats * Mats');          % [Nu_UE x Nu_UE], real

    %% ------------------------------------------------------------------
    % Frobenius norms and their outer-product matrix
    % ------------------------------------------------------------------
    norms   = sqrt(max(diag(GM), 0));     % [Nu_UE x 1]
    normMat = norms * norms';             % [Nu_UE x Nu_UE]

    %% ------------------------------------------------------------------
    % CMD
    % ------------------------------------------------------------------
    D_CMD = 1 - GM ./ (normMat + eps);

    %% ------------------------------------------------------------------
    % Enforce symmetry, zero diagonal, non-negativity
    % ------------------------------------------------------------------
    D_CMD = 0.5 * (D_CMD + D_CMD');
    D_CMD(1:Nu_UE+1:end) = 0;
    D_CMD = max(D_CMD, 0);

end

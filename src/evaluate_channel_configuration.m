function [CTn, TWn, ERR, gainVec, dvec_dB] = evaluate_channel_configuration( ...
        Htot, R_xy, supervised_Len, ArrSize, Nreal, dissimFn, ...
        tsne_perplexity, K_fraction)
% EVALUATE_CHANNEL_CONFIGURATION
% Computes channel-charting quality metrics for a single combined channel.
%
% This is the single-configuration counterpart of run_channel_charting,
% which loops over a full codebook.  Use this function to evaluate any
% individual channel hypothesis: a static EMS codeword, a dynamic/random
% baseline, the direct-path-only case, etc.
%
% Inputs:
%   Htot             : [Len x ArrSize x Nsc] combined channel (any config)
%   R_xy             : [Len x 2] ground-truth (x,y) user positions
%   supervised_Len   : number of labeled (supervised) users;
%                      they occupy rows 1:supervised_Len of Htot / R_xy
%   ArrSize          : number of BS antennas (= size(Htot,2))
%   Nreal            : number of i.i.d. noise realisations for covariance
%   dissimFn         : function handle for the dissimilarity metric,
%                      e.g.  @CompDismilarity_CMD
%                        or  @CompDismilarityLogEuclidean
%   tsne_perplexity  : (optional) t-SNE perplexity — default 150
%   K_fraction       : (optional) K = K_fraction * Len used for CT/TW
%                      — default 0.05
%
% Outputs:
%   CTn      : [(Len-supervised_Len) x 1]  −continuity  (unlabeled users)
%   TWn      : [(Len-supervised_Len) x 1]  −trustworthiness (unlabeled)
%   ERR      : [(Len-supervised_Len) x 1]  positioning errors [m]
%   gainVec  : [Len x 1]  channel gain per user [dB], 10*log10(trace(H*H'))
%   dvec_dB  : [N_pairs x 1]  upper-triangular pairwise dissimilarity [dB],
%              i.e. 10*log10(D_ij) for all i < j
%
% Convention note:
%   CTn and TWn are returned NEGATED (sign matches CTnegCell / TWnegCell
%   in the experiment mains) so that values closer to 0 indicate better
%   neighbourhood preservation.

    if nargin < 7 || isempty(tsne_perplexity), tsne_perplexity = 150;  end
    if nargin < 8 || isempty(K_fraction),      K_fraction      = 0.05; end

    Len = size(Htot, 1);
    KK  = K_fraction * Len;

    %% ------------------------------------------------------------------
    % 1. Sample covariance matrices + per-user channel gain
    % ------------------------------------------------------------------
    R       = zeros(Len, ArrSize, ArrSize);
    gainVec = zeros(Len, 1);

    for uu = 1:Len
        H = squeeze(Htot(uu, :, :));
        gainVec(uu) = 10 * log10(real(trace(H * H')) + eps);

        C = zeros(ArrSize);
        for r = 1:Nreal
            Hs = H * (randn + 1j*randn) / sqrt(2);
            C  = C + Hs * Hs';
        end
        R(uu, :, :) = C / Nreal;
    end

    %% ------------------------------------------------------------------
    % 2. Pairwise dissimilarity matrix
    % ------------------------------------------------------------------
    [D, ~]  = dissimFn(R);
    dvec_dB = 10 * log10(max(D(triu(true(size(D)), 1)), 1e-14));

    %% ------------------------------------------------------------------
    % 3. Semi-supervised t-SNE  (labeled users fixed to R_xy anchors)
    % ------------------------------------------------------------------
    xy = tsne_sem(D, R_xy, 1:supervised_Len, tsne_perplexity);

    %% ------------------------------------------------------------------
    % 4. Positioning error — unlabeled users only
    % ------------------------------------------------------------------
    ERR = sqrt(sum((xy - R_xy).^2, 2));
    ERR = ERR(supervised_Len+1:end);

    %% ------------------------------------------------------------------
    % 5. Continuity / Trustworthiness — unlabeled users only
    %    drscore_pointwise returns [Trustworthiness, Continuity] per row.
    %    We negate both to match the CTnegCell / TWnegCell convention
    %    used in the experiment mains (closer to 0 = better).
    % ------------------------------------------------------------------
    DistT = squareform(pdist(R_xy));
    DistZ = squareform(pdist(xy));
    TC    = drscore_pointwise(DistT, DistZ, KK, 'TC');

    CTn = -TC(supervised_Len+1:end, 1);   % −Trustworthiness (col 1)
    TWn = -TC(supervised_Len+1:end, 2);   % −Continuity      (col 2)

end

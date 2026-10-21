function F = abl_state_features(Data, params, st, noise)
% ABL_STATE_FEATURES  Log-covariance features of every user for one acquisition.
%   st.direct : include the direct BS-UE channel (H_LOS_1bounce)
%   st.phi    : [1x2] phase gradient of EMS1, EMS2 [deg]; NaN = surface off (ideal absorber)
%   st.key    : unique name of the acquisition (seeds its random numbers)
%   noise     : .name, .tx_dBm, .noise_dBm, .pilot_gain_dB (NaN = noiseless)
%   F         : [U x Nb^2] real, with ||F_u - F_v||_2 = ||log R_u - log R_v||_F (LED)
%
% Covariance per user, normalized by the UE transmit power sigma_s^2 (N = params.N_real snapshots,
% g_r ~ CN(0,1) as in the original code):
%   noiseless : R = (1/N) sum_r |g_r|^2 H H' + 2 sigma_n^2/sigma_s^2 I
%               (the floor only keeps R invertible; same level as the noise floor at 0 dB pilot gain)
%   estimated : R = (1/N) sum_r (H g_r + E_r)(H g_r + E_r)',  E_r ~ CN(0, sigma_e^2 I),
%               sigma_e^2 = sigma_n^2/sigma_s^2 / pilot gain.  Sampled exactly without the N loop:
%               projecting the E_r on g gives (||g|| H + F)(||g|| H + F)' + W, F ~ CN(0, sigma_e^2),
%               W ~ complex Wishart(2(N-1), sigma_e^2 I).

    nU = size(Data.R_xy,1);
    Nb = params.Par.Tx.ArrSize;
    Nr = size(Data.H_LOS_1bounce,3);

    if st.direct
        H = permute(double(Data.H_LOS_1bounce), [2 3 1]);          % Nb x Nr x U
    else
        H = complex(zeros(Nb, Nr, nU));
    end
    for e = 1:2
        if isnan(st.phi(e)), continue; end
        prof = kron(deg2rad(st.phi(e))*(1:params.M), ones(1,params.N));
        th   = exp(-1j*prof).';
        Ho   = double(Data.(sprintf('H_o_RIS%d',e)));
        Hi   = Data.(sprintf('H_I_RIS%d',e));                       % U x M*N x Nr
        X    = th .* reshape(permute(Hi,[2 3 1]), size(Hi,2), []);  % M*N x (Nr*U)
        H    = H + reshape(Ho*X, Nb, Nr, nU);
    end

    key = [st.key '|' noise.name];
    rs  = RandStream('mt19937ar', 'Seed', mod(sum(double(key).*(1:numel(key))), 2^31-2) + 1);
    sn2 = 10^((noise.noise_dBm - noise.tx_dBm)/10);
    noisy = ~isnan(noise.pilot_gain_dB);
    if noisy, se2 = sn2 / 10^(noise.pilot_gain_dB/10); end
    Nre = params.N_real;

    iu = find(triu(true(Nb),1));
    F  = zeros(nU, Nb + 2*numel(iu));
    for u = 1:nU
        Hu = H(:,:,u);
        v  = randn(rs, 2*Nre, 1);
        g2 = sum(v.^2)/2;                                  % sum_r |g_r|^2
        if noisy
            A = sqrt(g2)*Hu + sqrt(se2/2)*complex(randn(rs,Nb,Nr), randn(rs,Nb,Nr));
            R = (A*A' + se2*cwishart(rs, Nb, Nr*(Nre-1))) / Nre;
        else
            R = (g2/Nre)*(Hu*Hu') + 2*sn2*eye(Nb);
        end
        R = (R + R')/2;
        [V, ev] = eig(R, 'vector');
        L = (V .* log(ev).') * V';
        F(u,:) = [real(diag(L)).', sqrt(2)*real(L(iu)).', sqrt(2)*imag(L(iu)).'];
    end
end

function W = cwishart(rs, p, n)
% complex Wishart CW_p(n, I) by the Bartlett decomposition (E[W] = n I)
    T = zeros(p);
    for i = 1:p
        T(i,i) = sqrt(sum(randn(rs, 2*(n-i+1), 1).^2)/2);
    end
    lo = find(tril(true(p),-1));
    T  = complex(T);
    T(lo) = complex(randn(rs,numel(lo),1), randn(rs,numel(lo),1))/sqrt(2);
    W = T*T';
end

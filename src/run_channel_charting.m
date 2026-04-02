function results = run_channel_charting(state, params, codebook)

    Len = state.Len;
    nCodebooks = codebook.nCodebooks;

    %% PRE-ALLOCATION
    results.PosEstimates = zeros(Len,2,nCodebooks);
    results.ErrorPerCodebook = zeros(nCodebooks,1);
    results.ErrorPerUser = zeros(Len,nCodebooks);

    results.CT_PerUser = zeros(Len,nCodebooks);
    results.TW_PerUser = zeros(Len,nCodebooks);

    results.CT_PerCodebook_mean = zeros(nCodebooks,1);
    results.CT_PerCodebook_p90 = zeros(nCodebooks,1);
    results.TW_PerCodebook_mean = zeros(nCodebooks,1);
    results.TW_PerCodebook_p90 = zeros(nCodebooks,1);

    %% METRICS
    KK = params.K_fraction * Len;
    DistT = squareform(pdist(state.R_xy));

    fprintf('\n=== PHASE 1: Channel Charting for Each Codebook ===\n');

    %% LOOP
    for cb = 1:nCodebooks

        phi1 = codebook.all(cb,1);
        phi2 = codebook.all(cb,2);

        if mod(cb,10)==0 || cb==1
            fprintf('[Codebook %3d/%3d] φ1=%.4f rad, φ2=%.4f rad (%.1f°, %.1f°)\n', ...
                cb, nCodebooks, phi1, phi2, rad2deg(phi1), rad2deg(phi2));
        end

        %% ======= EXACT ORIGINAL LOGIC STARTS HERE =======

        % Phase profiles
        phase_profile1 = kron(phi1*(1:params.M), ones(1,params.N));
        phase_profile2 = kron(phi2*(1:params.M), ones(1,params.N));

        Theta_EMS1 = diag(exp(-1j*phase_profile1));
        Theta_EMS2 = diag(exp(-1j*phase_profile2));

        % Channel construction
        H_EMS1 = zeros(Len, params.Par.Tx.ArrSize, size(state.H_i_EMS1,3));
        H_EMS2 = zeros(Len, params.Par.Tx.ArrSize, size(state.H_i_EMS2,3));

        for uu = 1:Len
            % % % % disp(uu);
            H_EMS1(uu,:,:) = state.H_o_EMS1 * Theta_EMS1 * squeeze(state.H_i_EMS1(uu,:,:));
            H_EMS2(uu,:,:) = state.H_o_EMS2 * Theta_EMS2 * squeeze(state.H_i_EMS2(uu,:,:));
        end

        H_tot = state.H_dir + H_EMS1 + H_EMS2;

        % Covariance
        R_Cov = zeros(Len, params.Par.Tx.ArrSize, params.Par.Tx.ArrSize);

        for uu = 1:Len
            H = squeeze(H_tot(uu,:,:));
            C = zeros(params.Par.Tx.ArrSize);

            for r = 1:params.N_real
                Hs = H*(randn+1j*randn)/sqrt(2);
                C  = C + Hs*Hs';
            end

            R_Cov(uu,:,:) = C / params.N_real;
        end

        %% ======= EXACT ORIGINAL LOGIC ENDS HERE =======

        % Dissimilarity
        [D_mat, ~] = CompDismilarityLogEuclidean(R_Cov);

        % t-SNE
        xy_tsne = tsne_sem(D_mat, state.R_xy, ...
                           1:params.supervised_Len, ...
                           params.tsne_perplexity);

        results.PosEstimates(:,:,cb) = xy_tsne;

        % Error
        err = sqrt(sum((xy_tsne - state.R_xy).^2,2));
        results.ErrorPerUser(:,cb) = err;
        results.ErrorPerCodebook(cb) = mean(err(params.supervised_Len+1:end));

        % TW / CT
        DistZ = squareform(pdist(xy_tsne));
        TC = drscore_pointwise(DistT, DistZ, KK, 'TC');

        results.CT_PerUser(:,cb) = TC(:,1);
        results.TW_PerUser(:,cb) = TC(:,2);

        ct_unlabeled = TC(params.supervised_Len+1:end,1);
        tw_unlabeled = TC(params.supervised_Len+1:end,2);

        results.CT_PerCodebook_mean(cb) = mean(ct_unlabeled);
        results.CT_PerCodebook_p90(cb)  = prctile(ct_unlabeled,90);

        results.TW_PerCodebook_mean(cb) = mean(tw_unlabeled);
        results.TW_PerCodebook_p90(cb)  = prctile(tw_unlabeled,90);

        if mod(cb,10)==0 || cb==1
            fprintf('  → Avg error: %.2f m, Mean CT: %.3f, Mean TW: %.3f\n', ...
                results.ErrorPerCodebook(cb), ...
                results.CT_PerCodebook_mean(cb), ...
                results.TW_PerCodebook_mean(cb));
        end

    end

    fprintf('\nChannel charting completed for all %d codebooks.\n', nCodebooks);

end
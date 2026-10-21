function D = abl_job_input(job, getF)
% ABL_JOB_INPUT  Input dissimilarity matrix of one t-SNE job.
%   getF(direct, phi) returns the log-covariance features of all users for one acquisition.
%   F0 : features of the single configuration
%   F1 : configurations stacked
%   F2 : configurations stacked, each minus the user's mean over its centring group
%   F3 : configurations stacked, each minus the matching OFF acquisition
    u = job.chart.users;
    S = size(job.phi,1);
    X = [];
    for s = 1:S
        Fs = getF(job.direct, job.phi(s,:));
        X  = cat(3, X, Fs(u,:));
    end
    switch job.feature
        case 'F2'
            for g = unique(job.group).'
                m = job.group == g;
                X(:,:,m) = X(:,:,m) - mean(X(:,:,m), 3);
            end
        case 'F3'
            for s = 1:S
                Fo = getF(job.direct, job.off(s,:));
                X(:,:,s) = X(:,:,s) - Fo(u,:);
            end
    end
    X = reshape(X, numel(u), []);
    D = squareform(pdist(X));
end

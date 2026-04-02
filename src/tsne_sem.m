function ydata = tsne_sem(D, xy_true, selected_idx, perplexity)
%TSNE_P Performs symmetric t-SNE on affinity matrix P
%
%   mappedX = tsne_p(P, labels, no_dims)
%
% The function performs symmetric t-SNE on pairwise similarity matrix P 
% to create a low-dimensional map of no_dims dimensions (default = 2).
% The matrix P is assumed to be symmetric, sum up to 1, and have zeros
% on the diagonal.
% The labels of the data are not used by t-SNE itself, however, they 
% are used to color intermediate plots. Please provide an empty labels
% matrix [] if you don't want to plot results during the optimization.
% The low-dimensional data representation is returned in mappedX.
%
%

% First check whether we already have an initial solution
% set labels
    labels = xy_true(selected_idx,:);
    if numel(labels) > 1
        initial_solution = true;
        ydata_l   = labels;
        no_dims = size(ydata_l, 2);
        
    else
        initial_solution = false;
    end
    
    % Initialize some variables
    P = d2p(D, perplexity, 1e-5);  
    
    n = size(P, 1);                                     % number of instances
    momentum = 0.6;                                     % initial momentum .6
    final_momentum = 0.8;                               % value to which momentum is changed .8
    mom_switch_iter = 250;                              % iteration at which momentum is changed
    stop_lying_iter = 100;                              % iteration at which lying about P-values is stopped
    max_iter = 2000;                                    % maximum number of iterations
    epsilon  = 1000;                                     % initial learning rate
    min_gain = .01;                                     % minimum gain for delta-bar-delta
    
    % Make sure P-vals are set properly
    P(1:n + 1:end) = 0;                                 % set diagonal to zero
    P = 0.5 * (P + P');                                 % symmetrize P-values
    P = max(P ./ sum(P(:)), realmin);                   % make sure P-values sum to one
    const = sum(P(:) .* log(P(:)));                     % constant in KL divergence
    if initial_solution                                % ~initial_solution
        P = P * 4;                                      % lie about the P-vals to find better local minima
    end
    
    % Initialize the solution
    if ~initial_solution
        ydata = .0001 * randn(n, no_dims);              % used for unsupervised 
    else
        ydata = zeros(n, no_dims);                      % used for semi-supervised
        W_D = D;
        for i = 1:length(W_D)
            W_D(i,i) = inf;
            if ismember(i, selected_idx)
                ydata(i,:) = xy_true(i,:);
            else
                [val_x,idx] = min(W_D(i,selected_idx));
                ydata (i,:) = xy_true(selected_idx(idx),:);
            end
        end

    end
    
    y_incs  = zeros(size(ydata));
    gains   = ones(size(ydata));
    
    % Run the iterations solving KL func. using gradient descent
    for iter=1:max_iter
        
        % Compute joint probability that point i and j are neighbors
        sum_ydata = sum(ydata .^ 2, 2);
        num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata',...
            - 2 * (ydata * ydata')))); % Student-t distribution
        num(1:n+1:end) = 0;                                                 % set diagonal to zero
        Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
        
        % Compute the gradients (faster implementation)
        L = (P - Q) .* num;
        y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
            
        % Update the solution
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...         % note that the y_grads are actually -y_grads
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
        ydata = ydata + y_incs;
        % ydata = bsxfun(@minus, ydata, mean(ydata, 1));
        ydata(selected_idx,:) = labels;

        % Update the momentum if necessary
        if iter == mom_switch_iter
            momentum = final_momentum;
        end
        
        if iter == stop_lying_iter && initial_solution                     % ~initial_solution
            P = P ./ 4;
        end
        
        
        if ~rem(iter, 10)
            cost = const - sum(P(:) .* log(Q(:)));
            
            if iter>10 && abs(pre_cost - cost)/cost < 1e-7
                break;    %  stop!
            end
            %disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
            pre_cost = cost;
        end
        %%% moving points during CC
        % if ~rem(iter, 50)
        %     figure,
        %     scatter(ydata(:,1),ydata(:,2),50,ydata(:,2),'filled')
        %     hold on 
        %     scatter(labels(:,1),labels(:,2),50,'k')
        %     title(['number of iteration=',num2str(iter)])
        %     text(ydata(end-1:end,1), ydata(end-1:end,2), num2str(i), 'VerticalAlignment',...
        %         'bottom', 'HorizontalAlignment', 'right'); % 2 points are moving
        %     axis equal
        % end
        
       
        
        
    end
    % disp(iter)
    % disp(['Iteration ' num2str(iter) ': cost ' num2str(cost)]);
end  
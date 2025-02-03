% BIC method for multi-factor rank selection 

%%%%%%%%%%%%%%%%%%%%%%%%%
% input of bic_def_uni: %
%%%%%%%%%%%%%%%%%%%%%%%%%

% X: semi-symmetric tensors
% rank_max: the largest value of rank you want to try. For example, when p
% is large, it takes long time to try all different possible ranks. As
% we only care about the low rank structure of tensor factors Vk, we
% can set rank_max = 5.
% K: number of layers
% u0: initialization
% tol, max_iter: tolerance level and maximun iteration number


%%%%%%%%%%%%%%%%%%%%%%%%%%
% output of bic_def_uni: %
%%%%%%%%%%%%%%%%%%%%%%%%%%

% bic_rx: best rank selected by minimizing BIC. bic_rx are vectors of ranks
% for each Vk, so bic_rx is of dimension K

function [bic_rx] = bic_def_uni(X, rank_max, K, u0, spectral_init, tol, max_iter, deflation)

    bic_rx = ones(K, 1); % K is the number of layers, which is prespecified
    X_iter = cell(K+1, 1);
    
    k = 1;
    X_iter{1} = X;
    
    while k < K+1

        bic = ones(1, rank_max);

        if (k > 1) && spectral_init
            u0 = init_sst(X_iter{k});
        end

        for i = 1:rank_max
            bic(i) = bic_sst_uni(i, X_iter{k}, u0, tol, max_iter);
        end

        [~, bic_min] = min(bic(:));
        bic_rx(k) = bic_min; % bic estimation at kth layer

        % deflation to next layer
        [hat_u, hat_V, ~, hat_X] = sst_single(X_iter{k}, u0, bic_rx(k), tol, max_iter);
        
        if deflation == 0 % subtract deflation
            X_iter{k+1} = X_iter{k} - hat_X;
        elseif deflation == 1 % project deflation
            sz_X = size(X);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
        elseif deflation == 2 % project defaltion on joint factor
            X_iter{k+1} = X_iter{k} - hat_X;
            sz_X = size(X);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            X_iter{k+1} = ttm(X_iter{k+1}, {uk}, 3);
        elseif deflation == 3 % project deflation on individual factor
            X_iter{k+1} = X_iter{k} - hat_X;
            sz_X = size(X);
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            X_iter{k+1} = ttm(X_iter{k+1}, {Vk, Vk}, [1, 2]);
        end

        u0 = hat_u;
        k = k + 1;

    end

    bic_rx = bic_rx';

end
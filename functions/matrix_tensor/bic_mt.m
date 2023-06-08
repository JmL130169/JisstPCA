% BIC method for rank selection in matrix_tensor case 

% input of bic_mt:

% X, Y: semi-symmetric tensor and matrix
% rank_max: the largest value of rank you want to try. For example, when p
% is large, it takes long time to try all different possible ranks. As
% we only care about the low rank structure of tensor factors Vk, we
% can set rank_max = 5.
% K: number of layers
% u0, lambda, tol, max_iter: the input needed in Jimt_single
% lambda: vector of length K, each element is the relative weight of X
% compared with Y
% tol, max_iter: tolerance level and maximun iteration number
% method: if method = 1, we do not assume X and Y are of the same rank.
% When method = 2, we assume they have the same rank


% output of bic_def:

% bic_rx: best rank selected by minimizing BIC. bic_rx is vector of rank
% for each Vk, and bic_rx, bic_ry are of dimension K

function [bic_rx] = bic_mt(X, Y, rank_max, K, u0, lambda, tol, max_iter, deflation)

    bic_rx = ones(K, 1); % K is the number of layers, which is prespecified
    X_iter = cell(K+1, 1);
    Y_iter = cell(K+1, 1);
    
    k = 1;
    X_iter{1} = X;
    Y_iter{1} = Y;

    while k < K+1
        bic = ones(rank_max, 1);

        for i = 1:rank_max
            bic(i) = bic_msst(i, X_iter{k}, Y_iter{k}, u0, lambda(k), tol, max_iter);
        end
        [~, bic_min] = min(bic);
        bic_rx(k) = bic_min; % BIC estimation for the kth layer, X

        [hat_u, hat_V, hat_w, hat_X, hat_Y, ~, ~] = Jimt_single(X_iter{k}, Y_iter{k}, u0, bic_rx(k), lambda(k), tol, max_iter);

        if deflation == 0 % subtract deflation
            X_iter{k+1} = X_iter{k} - hat_X;
            Y_iter{k+1} = Y_iter{k} - hat_Y;
        elseif deflation == 1 % project deflation
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            wk = double(eye(sz_Y(1)) - hat_w*hat_w');
            X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_iter{k+1} = wk*Y_iter{k}*uk;
        elseif deflation == 2 % project defaltion on joint factor
            X_iter{k+1} = X_iter{k} - hat_X;
            Y_iter{k+1} = Y_iter{k} - hat_Y;
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            wk = double(eye(sz_Y(1)) - hat_w*hat_w');
            X_iter{k+1} = ttm(X_iter{k+1}, {uk}, 3);
            Y_iter{k+1} = Y_iter{k+1}*uk;
        elseif deflation == 3 % project deflation on individual factor
            X_iter{k+1} = X_iter{k} - hat_X;
            Y_iter{k+1} = Y_iter{k} - hat_Y;
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            wk = double(eye(sz_Y(1)) - hat_w*hat_w');
            X_iter{k+1} = ttm(X_iter{k+1}, {Vk, Vk}, [1, 2]);
            Y_iter{k+1} = wk*Y_iter{k+1};
        end

        u0 = hat_u;

        k = k+1;
    end
    

    bic_rx = bic_rx';
    
end

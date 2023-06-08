% BIC method for multi-factor rank selection 

% input of bic_def:

% X, Y: semi-symmetric tensors
% rank_max: the largest value of rank you want to try. For example, when p,
% q are large, it takes long time to try all different possible ranks. As
% we only care about the low rank structure of tensor factors Vk and Wk, we
% can set rank_max = 5.
% K: number of layers
% u0, lambda, tol, max_iter: the input needed in Jisst_single
% lambda: vector of length K, each element is the relative weight of X
% compared with Y
% tol, max_iter: tolerance level and maximun iteration number
% method: if method = 1, we do not assume X and Y are of the same rank.
% When method = 2, we assume they have the same rank


% output of bic_def:

% bic_rx, bic_ry: best rank selected by minimizing BIC. bic_rx, bic_ry are vectors of ranks
% for each Vk, Wk, so bic_rx, bic_ry are of dimension K

function [bic_rx, bic_ry] = bic_def(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, deflation)

    bic_rx = ones(K, 1); % K is the number of layers, which is prespecified
    bic_ry = ones(K, 1);
    X_iter = cell(K+1, 1);
    Y_iter = cell(K+1, 1);
    
    k = 1;
    X_iter{1} = X;
    Y_iter{1} = Y;

    if method == 1
        while k < K+1
            bic = ones(rank_max, rank_max);

            for i = 1:rank_max
                for j = 1:rank_max
                    bic(i, j) = bic_sst(i, j, X_iter{k}, Y_iter{k}, u0, lambda(k), tol, max_iter);
                end
            end
            [~, bic_min] = min(bic(:));
            [bic_rx(k), bic_ry(k)] = ind2sub(size(bic), bic_min); % bic estimation at kth layer

            [hat_u, hat_V, hat_W, ~, ~, hat_X, hat_Y] = Jisst_single(X_iter{k}, Y_iter{k}, u0, bic_rx(k), bic_ry(k), lambda(k), tol, max_iter);
            
            if deflation == 0 % subtract deflation
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
            elseif deflation == 1 % project deflation
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
                Y_iter{k+1} = ttm(Y_iter{k}, {Wk, Wk, uk}, [1, 2, 3]);
            elseif deflation == 2 % project defaltion on joint factor
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k+1}, {uk}, 3);
                Y_iter{k+1} = ttm(Y_iter{k+1}, {uk}, 3);
            elseif deflation == 3 % project deflation on individual factor
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k+1}, {Vk, Vk}, [1, 2]);
                Y_iter{k+1} = ttm(Y_iter{k+1}, {Wk, Wk}, [1, 2]);
            end

            u0 = hat_u;

            k = k + 1;
        end

    elseif method == 2
        while k < K+1
            bic = ones(rank_max, 1);

            for i = 1:rank_max
                bic(i) = bic_sst(i, i, X_iter{k}, Y_iter{k}, u0, lambda(k), tol, max_iter);
            end
            [~, bic_min] = min(bic);
            bic_rx(k) = bic_min; % BIC estimation for the kth layer, X
            bic_ry(k) = bic_min; % BIC estimation for the kth layer, Y

            [hat_u, hat_V, hat_W, ~, ~, hat_X, hat_Y] = Jisst_single(X_iter{k}, Y_iter{k}, u0, bic_rx(k), bic_ry(k), lambda(k), tol, max_iter);
            
            if deflation == 0 % subtract deflation
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
            elseif deflation == 1 % project deflation
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
                Y_iter{k+1} = ttm(Y_iter{k}, {Wk, Wk, uk}, [1, 2, 3]);
            elseif deflation == 2 % project defaltion on joint factor
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k+1}, {uk}, 3);
                Y_iter{k+1} = ttm(Y_iter{k+1}, {uk}, 3);
            elseif deflation == 3 % project deflation on individual factor
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k+1}, {Vk, Vk}, [1, 2]);
                Y_iter{k+1} = ttm(Y_iter{k+1}, {Wk, Wk}, [1, 2]);
            end
            
            u0 = hat_u;

            k = k+1;
        end
    end

    bic_rx = bic_rx';
    bic_ry = bic_ry';
end


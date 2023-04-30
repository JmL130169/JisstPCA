% BIC method for multi-factor rank selection based on single-factor BIC
% with subtract deflation

% input of bic_def:
% X, Y: multi-factor semi-symmetric tensors
% rank_max: the largest value of rank you want to try. For example, when p,
% q are large, it takes long time to try all different possible ranks. As
% we only care about the low rank structure of tensor factors Vk and Wk, we
% can set rank_max = 5.
% K: number of layers
% u0, lambda, tol, max_iter: the input needed in Jisst_single
% method: if method = 1, we do not assume X and Y are of the same rank.
% When method = 2, we assume they have the same rank
% def: deflation strategy. def = 0 is subtract deflation; def = 1 is project deflation; def = 2 is subtract deflation for X and project deflation for Y; def = 3 is project deflation for X and subtract deflation for Y

% output of bic_def:
% bic_rx, bic_ry: best rank selected by minimizing BIC. bic_rx, bic_ry are vectors of ranks
% for each Vk, Wk, so bic_rx, bic_ry are of dimension K

% Note: this function takes long time. When applying this function, you
% might assume the rank of X equals the rank of Y by inputting method = 2

function [bic_rx, bic_ry] = bic_def(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, def)

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
                    bic(i, j) = bic_sst(i, j, X_iter{k}, Y_iter{k}, u0, lambda, tol, max_iter);
                end
            end
            [~, bic_min] = min(bic);
            bic_rx(k) = bic_min(1); % BIC estimation for the kth layer, X
            bic_ry(k) = bic_min(2); % BIC estimation for the kth layer, Y

            [hat_u, hat_V, hat_W, ~, ~, hat_X, hat_Y] = Jisst_single(X_iter{k}, Y_iter{k}, u0, bic_rx(k), bic_ry(k), lambda, tol, max_iter);
            
            if def == 0 % subtract deflation
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
            else if def == 1 % project deflation
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
                Y_iter{k+1} = ttm(Y_iter{k}, {Wk, Wk, uk}, [1, 2, 3]);
            else if def == 2 % subtract deflation for X and project deflation for Y
                X_iter{k+1} = X_iter{k} - hat_X;
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                Y_iter{k+1} = ttm(Y_iter{k}, {Wk, Wk, uk}, [1, 2, 3]);
            else if def == 3 % project deflation for X and subtract deflation for Y
                sz_X = size(X);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
                Y_iter{k+1} = Y_iter{k} - hat_Y;
            end
            
            u0 = hat_u;

            k = k + 1;
        end

    elseif method == 2
        while k < K+1
            bic = ones(rank_max, 1);

            for i = 1:rank_max
                bic(i) = bic_sst(i, i, X_iter{k}, Y_iter{k}, u0, lambda, tol, max_iter);
            end
            [~, bic_min] = min(bic);
            bic_rx(k) = bic_min; % BIC estimation for the kth layer, X
            bic_ry(k) = bic_min; % BIC estimation for the kth layer, Y

            [hat_u, hat_V, hat_W, ~, ~, hat_X, hat_Y] = Jisst_single(X_iter{k}, Y_iter{k}, u0, bic_rx(k), bic_ry(k), lambda, tol, max_iter);
            
            if def == 0 % subtract deflation
                X_iter{k+1} = X_iter{k} - hat_X;
                Y_iter{k+1} = Y_iter{k} - hat_Y;
            else if def == 1 % project deflation
                sz_X = size(X);
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
                Y_iter{k+1} = ttm(Y_iter{k}, {Wk, Wk, uk}, [1, 2, 3]);
            else if def == 2 % subtract deflation for X and project deflation for Y
                X_iter{k+1} = X_iter{k} - hat_X;
                sz_Y = size(Y);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Wk = double(eye(sz_Y(1)) - hat_W*hat_W');
                Y_iter{k+1} = ttm(Y_iter{k}, {Wk, Wk, uk}, [1, 2, 3]);
            else if def == 3 % project deflation for X and subtract deflation for Y
                sz_X = size(X);
                uk = double(eye(sz_X(3)) - hat_u*hat_u');
                Vk = double(eye(sz_X(1)) - hat_V*hat_V');
                X_iter{k+1} = ttm(X_iter{k}, {Vk, Vk, uk}, [1, 2, 3]);
                Y_iter{k+1} = Y_iter{k} - hat_Y;
            end
            
            u0 = hat_u;

            k = k+1;
        end
    end
end

% Note: the output of bic_def might not be directly used in Jisst_multi or
% dJisst_multi because bic_rx, bic_ry we obtained are column vectors. We
% might need to simply transpose the output before using them.

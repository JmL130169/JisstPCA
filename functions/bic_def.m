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

% output of bic_def:
% bic_rx, bic_ry: best rank selected by minimizing BIC. bic_rx, bic_ry are vectors of ranks
% for each Vk, Wk, so bic_rx, bic_ry are of dimension K

% Note: this function takes long time. When applying this function, you
% might assume the rank of X equals the rank of Y by inputting method = 2

function [bic_rx, bic_ry] = bic_def(X, Y, rank_max, K, u0, lambda, tol, max_iter, method)

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
            [~, bic_min] = min(bic(:));
            [bic_rx(k), bic_ry(k)] = ind2sub(size(bic), bic_min); % bic estimation at kth layer

            rank{1} = bic_rx(k); rank{2} = bic_ry(k);
            r_max = 0;
            [hat_u, ~, ~, ~, ~, hat_X, hat_Y] = Jisst_single(X_iter{k}, Y_iter{k}, u0, lambda, tol, max_iter, r_max, rank);
            X_iter{k+1} = X_iter{k} - hat_X;
            Y_iter{k+1} = Y_iter{k} - hat_Y; % subtract deflation of the next layer
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

            [hat_u, ~, ~, ~, ~, hat_X, hat_Y] = Jisst_single(X_iter{k}, Y_iter{k}, u0, bic_rx(k), bic_ry(k), lambda, tol, max_iter);
            X_iter{k+1} = X_iter{k} - hat_X;
            Y_iter{k+1} = Y_iter{k} - hat_Y; % subtract deflation of the next layer
            u0 = hat_u;

            k = k+1;
        end
    end
end

% Note: the output of bic_def might not be directly used in Jisst_multi or
% dJisst_multi because bic_rx, bic_ry we obtained are column vectors. We
% might need to simply transpose the output before using them.

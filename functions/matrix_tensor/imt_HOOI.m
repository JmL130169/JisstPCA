% integrated Matrix-Tensor High-Order Orthogonal Iteration (iMT-HOOI)
% This function is similar to tensor-tensor version iHOOI

% input of imt_HOOI:
% X is semi-symmetric tensor with dimension p-p-N
% Y is rank K matrix of dimension q-N
% rx is the sum of rank of all factors V_k for X
% max_iter is the maximum iteration number
% K is the number of layers

% The output of imt_HOOI:
% hat_u is N-K matrix, ith column is the estimation of u_i
% hat_V is p-rx matrix, which can be viewed estimation of concatenated
% factors V_k, as hat_V = (hat_V1, hat_V2, ..., hat_VK)
% hat_W is q-K matrix, ith column is the estimation of w_i
% hat_X, hat_Y are the corresponding estimators of X and Y

function [hat_u, hat_V, hat_w, hat_X, hat_Y] = imt_HOOI(X, Y, rx, max_iter, K)

    % rx is sum of rank for each layer, which is not vector

    % initialization
    M1_X = double(tenmat(X, 1, [2, 3]));
    [ax, ~, ~] = svds(M1_X, rx);
    hat_V = ax;

    [ay, ~, ~] = svds(Y, K);
    hat_w = ay;

    M3_X = double(tenmat(X, 3, [1, 2]));
    M3 = [M3_X, Y'];
    [a, ~, ~] = svds(M3, K);
    hat_u = a;
    k = 0;

    % update of V, W, u
    while k < max_iter

        T_x = ttm(X, {hat_V', hat_u'}, [2, 3]);

        M1_X = double(tenmat(T_x, 1, [2, 3]));
        [ax, ~, ~] = svds(M1_X, rx);
        hat_V = ax;

        T_y = Y*hat_u;

        [ay, ~, ~] = svds(T_y, K);
        hat_w = ay;


        T_xu = ttm(X, {hat_V', hat_V'}, [1, 2]);
        M3_X = double(tenmat(T_xu, 3, [1, 2]));
        T_yu = hat_w'*Y;

        M3 = [M3_X, T_yu'];
        [a, ~, ~] = svds(M3, K);
        hat_u = a;

        k = k+1;
    end

    % reconstruct the true tensor X and Y
    S_x = ttm(X, {hat_V', hat_V', hat_u'}, [1, 2, 3]);
    D_y = hat_w'*Y*hat_u;
    hat_X = ttm(S_x, {hat_V, hat_V, hat_u}, [1, 2, 3]);
    hat_Y = hat_w*D_y*hat_u';
end

% Note: finally, the estimation of Vk is hat_V(:, rx_cum(k-1)+1:rx_cum(k)),
% where rx_cum is the cumsum of the vector of ranks of X.



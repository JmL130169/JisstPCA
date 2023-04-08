% integrated Matrix-Tensor High-Order SVD (iMT-HOSVD)
% This function is similar to tensor-tensor version iHOSVD

% input of imt_HOSVD:
% X is semi-symmetric tensor with dimension p-p-N
% Y is rank K matrix of dimension q-N
% rx is the sum of rank of all factors V_k for X
% nl = K is the number of layers

% The output of imt_HOSVD:
% hat_u is N-K matrix, ith column is the estimation of u_i
% hat_V is p-rx matrix, which can be viewed estimation of concatenated
% factors V_k, as hat_V = (hat_V1, hat_V2, ..., hat_VK)
% hat_W is q-K matrix, ith column is the estimation of w_i
% hat_X, hat_Y are the corresponding estimators of X and Y

function [hat_u, hat_V, hat_w, hat_X, hat_Y] = imt_HOSVD(X, Y, rx, nl)

    % nl is the number of layers, which is K in the model
    % rx is sum of rank for each layer, which is not vector
    M1_X = double(tenmat(X, 1, [2, 3]));
    M3 = [double(tenmat(X, 3, [1, 2])), Y'];

    [ax, bx, ~] = svd(M1_X);
    [~, ind_x] = sort(diag(bx));
    ax = ax(:, ind_x);
    ax = fliplr(ax);
    hat_V = ax(:, 1:rx); % estimation of V in Tucker model

    [ay, by, ~] = svd(Y);
    [~, ind_y] = sort(diag(by));
    ay = ay(:, ind_y);
    ay = fliplr(ay);
    hat_w = ay(:, 1:nl); % estimation of w in Tucker model

    [a, b, ~] = svd(M3);
    [~, ind] = sort(diag(b));
    a = a(:, ind);
    a = fliplr(a);
    hat_u = a(:, 1:nl); % estimation of u in Tucker model

    % reconstruct the true tensor X and Y
    S_x = ttm(X, {hat_V', hat_V', hat_u'}, [1, 2, 3]);
    D_y = hat_w'*Y*hat_u;
    hat_X = ttm(S_x, {hat_V, hat_V, hat_u}, [1, 2, 3]);
    hat_Y = hat_w*D_y*hat_u';
end

% Note: finally, the estimation of Vk is hat_V(:, rx_cum(k-1)+1:rx_cum(k)),
% where rx_cum is the cumsum of the vector of ranks of X.
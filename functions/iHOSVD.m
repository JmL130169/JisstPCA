% This function is high-order SVD for integrated
% semi-symmetric tensors (iHOSVD), by concatenated X and Y and then estimate
% each factor using classical HOSVD algorithm. Finally we put the output of
% iHOSVD to the corresponding factors in JisstPCA model by the relation
% between JisstPCA model and Tucker model

% The input of iHOSVD:
% X, Y are semi-symmetric tensors with dimension p-p-N and q-q-N, and
% suppose X, Y are multi-factor tensors with number of layers K
% rx is the sum of rank of all factors V_k for X, ry is the sum of rank of
% all factors W_k for Y
% nl = K is the number of layers

% The output of iHOSVD:
% hat_u is N-K matrix, ith column is the estimation of u_i
% hat_V is p-rx matrix, which can be viewed estimation of concatenated
% factors V_k, as hat_V = (hat_V1, hat_V2, ..., hat_VK)
% hat_W is q-ry matrix, which can be viewed estimation of concatenated
% factors W_k, as hat_W = (hat_W1, hat_W2, ..., hat_WK)
% hat_X, hat_Y are the corresponding estimators of X and Y

function [hat_u, hat_V, hat_W, hat_X, hat_Y] = iHOSVD(X, Y, rx, ry, nl)
    
    M1_X = double(tenmat(X, 1, [2, 3]));
    M1_Y = double(tenmat(Y, 1, [2, 3]));
    M3 = [double(tenmat(X, 3, [1, 2])), double(tenmat(Y, 3, [1, 2]))];

    % estimation of V is first rx singular vectors of first-mode
    % matricization of X
    [ax, bx, ~] = svd(M1_X);
    [~, ind_x] = sort(diag(bx));
    ax = ax(:, ind_x);
    ax = fliplr(ax);
    hat_V = ax(:, 1:rx);

    % estimation of W is first ry singular vectors of first-mode
    % matricization of Y
    [ay, by, ~] = svd(M1_Y);
    [~, ind_y] = sort(diag(by));
    ay = ay(:, ind_y);
    ay = fliplr(ay);
    hat_W = ay(:, 1:ry);
    
    % estimation of u is first K singular vectors of concatenated
    % matrix of third-mode matricization of X and third-mode matricization
    % of Y
    [a, b, ~] = svd(M3);
    [~, ind] = sort(diag(b));
    a = a(:, ind);
    a = fliplr(a);
    hat_u = a(:, 1:nl);

    % reconstruct the true tensor X and Y
    S_x = ttm(X, {hat_V', hat_V', hat_u'}, [1, 2, 3]);
    S_y = ttm(Y, {hat_W', hat_W', hat_u'}, [1, 2, 3]);
    hat_X = ttm(S_x, {hat_V, hat_V, hat_u}, [1, 2, 3]);
    hat_Y = ttm(S_y, {hat_W, hat_W, hat_u}, [1, 2, 3]);
end

% Note: finally, the estimation of Vk is hat_V(:, rx_cum(k-1)+1:rx_cum(k)),
% where rx_cum is the cumsum of the vector of ranks of X. And so is for Y
% and hat_W, Wk
% This function is high-order orthogonal iteration for integrated
% semi-symmetric tensors (iHOOI), by concatenated X and Y and then estimate
% each factor using classical HOOI algorithm. Finally we put the output of
% iHOOI to the corresponding factors in JisstPCA model by the relation
% between JisstPCA model and Tucker model

% The input of iHOOI:

% X, Y are semi-symmetric tensors with dimension p-p-N and q-q-N, and
% suppose X, Y are multi-factor tensors with number of layers K
% rx is the sum of rank of all factors V_k for X, ry is the sum of rank of
% all factors W_k for Y
% max_iter is the maximum iteration number
% K is the number of layers

% The output of iHOOI:

% hat_u is N-K matrix, ith column is the estimation of u_i
% hat_V is p-rx matrix, which can be viewed estimation of concatenated
% factors V_k, as hat_V = (hat_V1, hat_V2, ..., hat_VK)
% hat_W is q-ry matrix, which can be viewed estimation of concatenated
% factors W_k, as hat_W = (hat_W1, hat_W2, ..., hat_WK)
% hat_X, hat_Y are the corresponding estimators of X and Y

function [hat_u, hat_V, hat_W, hat_X, hat_Y] = iHOOI(X, Y, rx, ry, max_iter, K)

    % initialization

    % initialization of V is first rx singular vectors of first-mode
    % matricization of X
    M1_X = double(tenmat(X, 1, [2, 3]));
    [ax, ~, ~] = svds(M1_X, rx);
    hat_V = ax;
    
    % initialization of W is first ry singular vectors of first-mode
    % matricization of Y
    M1_Y = double(tenmat(Y, 1, [2, 3]));
    [ay, ~, ~] = svds(M1_Y, ry);
    hat_W = ay;
    
    % initialization of u is first K singular vectors of concatenated
    % matrix of third-mode matricization of X and third-mode matricization
    % of Y
    M3_X = double(tenmat(X, 3, [1, 2]));
    M3_Y = double(tenmat(Y, 3, [1, 2]));
    M3 = [M3_X, M3_Y];
    [a, ~, ~] = svds(M3, K);
    hat_u = a;

    k = 0;

    % update of V, W, u
    while k < max_iter

        T_x = ttm(X, {hat_V', hat_u'}, [2, 3]);

        M1_X = double(tenmat(T_x, 1, [2, 3]));
        [ax, ~, ~] = svds(M1_X, rx);
        hat_V = ax;


        T_y = ttm(Y, {hat_W', hat_u'}, [2, 3]);

        M1_Y = double(tenmat(T_y, 1, [2, 3]));
        [ay, ~, ~] = svds(M1_Y, ry);
        hat_W = ay;


        T_xu = ttm(X, {hat_V', hat_V'}, [1, 2]);
        M3_X = double(tenmat(T_xu, 3, [1, 2]));
        T_yu = ttm(Y, {hat_W', hat_W'}, [1, 2]);
        M3_Y = double(tenmat(T_yu, 3, [1, 2]));

        M3 = [M3_X, M3_Y];
        [a, ~, ~] = svds(M3, K);
        hat_u = a;

        k = k+1;
    end

    % reconstruct the true tensor X and Y
    S_x = ttm(X, {hat_V', hat_V', hat_u'}, [1, 2, 3]);
    S_y = ttm(Y, {hat_W', hat_W', hat_u'}, [1, 2, 3]);
    hat_X = ttm(S_x, {hat_V, hat_V, hat_u}, [1, 2, 3]);
    hat_Y = ttm(S_y, {hat_W, hat_W, hat_u}, [1, 2, 3]);
end

% Note: finally, the estimation of Vk is hat_V(:, rx_cum(k-1)+1:rx_cum(k)),
% where rx_cum is the cumsum of the vector of ranks of X. And so is for Y
% and hat_W, Wk


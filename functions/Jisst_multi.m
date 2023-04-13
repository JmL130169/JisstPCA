% This function is about multi-factor JisstPCA with subtract deflation

% Input of Jisst_multi is:
% X, Y: two multi-factor semi-symmetric tensors of dimension p-p-N and q-q-N, with K layers
% u0: initialization
% rx, ry: rank of X and Y, both are vectors of dimension K
% lambda: scaler for each layer, which is a vector of dimension K
% tol, max_iter: tolerance value and maximum iteration number

% Output of Jisst_multi is:
% u_est, V_est, W_est: cells with (K+1) elements. The estimation of kth factor of u, V, W are the (k+1)th factor in u_est, V_est, W_est
% d_est: matrix of dimension 2-(K+1), while the first row is estimation of dx, and second row is estimation of dy. The estimation of kth factor is in the (k+1)th column of d_est

function [u_est, V_est, W_est, d_est] = Jisst_multi(X, Y, u0, rx, ry, lambda, tol, max_iter)

    sz = size(rx);
    X_est = cell(sz(2) + 1, 1);
    Y_est = cell(sz(2) + 1, 1);
    u_est = cell(sz(2) + 1, 1);
    V_est = cell(sz(2) + 1, 1);
    W_est = cell(sz(2) + 1, 1);
    d_est = zeros(2, sz(2)+1);

    u_est{1} = u0;
    X_est{1} = X;
    Y_est{1} = Y;

    k = 1;
    while k < sz(2) + 1
        [hat_u, hat_V, hat_W, d_x, d_y, ~, ~] = Jisst_single(X_est{k}, Y_est{k}, u_est{k}, rx(k), ry(k), lambda(k), tol, max_iter);
        
        % update of tensor factors
        u_est{k+1} = hat_u;
        V_est{k+1} = hat_V;
        W_est{k+1} = hat_W;
        d_est(1, k+1) = d_x; % estimation of d_{x_{k}}
        d_est(2, k+1) = d_y; % estimation of d_{y_{k}}

        % deflate
        X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
        Y_est{k+1} = Y_est{k} - d_y*squeeze(ttt(tensor(hat_W*hat_W'), tensor(hat_u)));
        
        % finish the while loop
        k = k+1;
    end

    % by this function, the estimation of kth factor is (k+1)th element in
    % each cell, i.e. \hat{u}_{k} = u_est{k+1} and also for V, W
end

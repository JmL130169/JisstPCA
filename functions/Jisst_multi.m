% This function is about multi-factor JisstPCA with subtract deflation
% scheme
function [u_est, V_est, W_est, d_est] = Jisst_multi(X, Y, u0, rx, ry, lambda, tol, max_iter)

    % Input X and Y are multi-factor semi-symmetric tensors, rx and ry are
    % ranks of dimension K, and lambda is also a vector of dimension K
    % representing the scaler of each layer
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

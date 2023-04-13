% function of multi-factor JISST-PCA while the true model uses diagonal
% matrix D as signal, based on subtract deflation

% Input of dJisst_multi is:
% X, Y: two multi-factor semi-symmetric tensors of dimension p-p-N and q-q-N, with K layers
% u0: initialization
% rx, ry: rank of X and Y, both are vectors of dimension K
% lambda: scaler for each layer, which is a vector of dimension K
% tol, max_iter: tolerance value and maximum iteration number

% Output of dJisst_multi is:
% u_est, V_est, W_est: cells with (K+1) elements. The estimation of kth factor of u, V, W are the (k+1)th factor in u_est, V_est, W_est
% Dx_est, Dy_est: cells with (K+1) elements. The estimation of kth diagonal signal of Dx, Dy are the (k+1)th factor in Dx_est, Dy_est

function [u_est, V_est, W_est, Dx_est, Dy_est] = dJisst_multi(X, Y, u0, rx, ry, lambda, tol, max_iter)

    sz = size(rx, 2);
    X_est = cell(sz + 1, 1);
    Y_est = cell(sz + 1, 1);
    u_est = cell(sz + 1, 1);
    V_est = cell(sz + 1, 1);
    W_est = cell(sz + 1, 1);
    Dx_est = cell(sz + 1, 1);
    Dy_est = cell(sz + 1, 1);

    % set initialization and let X^{1} = X, Y^{1} = Y as in our algorithm
    u_est{1} = u0;
    X_est{1} = X;
    Y_est{1} = Y;

    k = 1;
    while k < (sz + 1)
        [hat_u, hat_V, hat_W, Dx, Dy, ~, ~] = dJisst_single(X_est{k}, Y_est{k}, u_est{k}, rx(k), ry(k), lambda(k), tol, max_iter);
        
        % update of tensor factors
        u_est{k+1} = hat_u;
        V_est{k+1} = hat_V;
        W_est{k+1} = hat_W;
        Dx_est{k+1} = Dx;
        Dy_est{k+1} = Dy;

        % subtract deflation
        X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
        Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));

        % next iteration
        k = k+1;
    end

    % Note: Again, to avoid confusion, the estimation of kth factor is
    % (k+1)th element of the output estimators, e.g. \hat{u}_{k} = u_est{k+1}
end








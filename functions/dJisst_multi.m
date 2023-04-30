% function of multi-factor JISST-PCA while the true model uses diagonal
% matrix D as signal, based on subtract deflation

% Input of dJisst_multi is:
% X, Y: two multi-factor semi-symmetric tensors of dimension p-p-N and q-q-N, with K layers
% u0: initialization
% rx, ry: rank of X and Y, both are vectors of dimension K
% lambda: scaler for each layer, which is a vector of dimension K
% tol, max_iter: tolerance value and maximum iteration number
% def: deflation strategy. def = 0 is subtract deflation; def = 1 is project deflation; def = 2 is subtract deflation for X and project deflation for Y; def = 3 is project deflation for X and subtract deflation for Y

% Output of dJisst_multi is:
% u_est, V_est, W_est: cells with (K+1) elements. The estimation of kth factor of u, V, W are the (k+1)th factor in u_est, V_est, W_est
% Dx_est, Dy_est: cells with (K+1) elements. The estimation of kth diagonal signal of Dx, Dy are the (k+1)th factor in Dx_est, Dy_est

function [u_est, V_est, W_est, Dx_est, Dy_est] = dJisst_multi(X, Y, u0, rx, ry, lambda, tol, max_iter, def)

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

        % deflation
        if def == 0 % subtract deflation
            X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));
        else if def == 1 % project deflation
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = ttm(Y_est{k}, {Wk, Wk, uk}, [1, 2, 3]);
        else if def == 2 % subtract deflation for X and project deflation for Y
            X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');           
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            Y_est{k+1} = ttm(Y_est{k}, {Wk, Wk, uk}, [1, 2, 3]);
        else if def == 3 % project deflation for X and subtract deflation for Y
            sz_X = size(X);            
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');           
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));
        end
        
        % next iteration
        k = k+1;
    end

    % Note: Again, to avoid confusion, the estimation of kth factor is
    % (k+1)th element of the output estimators, e.g. \hat{u}_{k} = u_est{k+1}
end








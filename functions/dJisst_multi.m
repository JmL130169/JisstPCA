% function of multi-factor JISST-PCA while the true model uses diagonal
% matrix D as signal, based on subtract deflation

% Input of dJisst_multi is:
% X, Y: two multi-factor semi-symmetric tensors of dimension p-p-N and q-q-N, with K layers
% u0: initialization
% rx, ry: rank of X and Y, both are vectors of dimension K
% lambda: scaler for each layer, which is a vector of dimension K
% tol, max_iter: tolerance value and maximum iteration number
% rank_max: if the true rank rx, ry are unknown, rank_max is the largest
% possible that will be tried using BIC method; If rx, ry are known, this
% arguemtn will not affect this function and user can set rank_max = 0
% def: def = 0 is subtract deflation, def = 1 is project deflation, def = 2
% is project deflation for only u after subtract deflation, def = 3 is
% project deflation for only V and W after subtract deflation
% varargin: If rx, ry are unknown, this should be a cell for the true
% ranks, if true ranks are unknown, this arguement can be left empty and
% JisstPCA will use BIC estimated ranks

% Output of dJisst_multi is:
% u_est, V_est, W_est: cells with (K+1) elements. The estimation of kth factor of u, V, W are the (k+1)th factor in u_est, V_est, W_est
% Dx_est, Dy_est: cells with (K+1) elements. The estimation of kth diagonal signal of Dx, Dy are the (k+1)th factor in Dx_est, Dy_est

function [u_est, V_est, W_est, Dx_est, Dy_est] = dJisst_multi(X, Y, u0, lambda, tol, max_iter, rank_max, def, varargin)

    % estimate rank if true rank is unknown
    if isempty(varargin)
        K = size(lambda, 2);
        method = 1;
        [bic_rx, bic_ry] = bic_def(X, Y, rank_max, K, u0, lambda, tol, max_iter, method);
        rx = bic_rx';
        ry = bic_ry';
    else
        r = varargin{1};
        rx = r{1};
        ry = r{2};
    end

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
        rank{1} = rx(k); rank{2} = ry(k);
        r_max = 0;
        [hat_u, hat_V, hat_W, Dx, Dy, ~, ~] = dJisst_single(X_est{k}, Y_est{k}, u_est{k}, lambda(k), tol, max_iter, r_max, rank);
        
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
        elseif def == 1 % project deflation
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = ttm(Y_est{k}, {Wk, Wk, uk}, [1, 2, 3]);
        elseif def == 2 % orthogonal joint factor after subtract deflation
            X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));
            
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');

            X_est{k+1} = ttm(X_est{k+1}, {uk}, 3);
            Y_est{k+1} = ttm(Y_est{k+1}, {uk}, 3);
        elseif def == 3 % orthogonal individual factors after subtract deflation
            X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));
            
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            
            X_est{k+1} = ttm(X_est{k+1}, {Vk, Vk}, [1, 2]);
            Y_est{k+1} = ttm(Y_est{k+1}, {Wk, Wk}, [1, 2]);
        end
        
        % next iteration
        k = k+1;
    end

    % Note: Again, to avoid confusion, the estimation of kth factor is
    % (k+1)th element of the output estimators, e.g. \hat{u}_{k} = u_est{k+1}
end








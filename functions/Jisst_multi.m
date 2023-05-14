% This function is about multi-factor JisstPCA with subtract deflation

% Input of Jisst_multi is:
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

% Output of Jisst_multi is:
% u_est, V_est, W_est: cells with (K+1) elements. The estimation of kth factor of u, V, W are the (k+1)th factor in u_est, V_est, W_est
% d_est: matrix of dimension 2-(K+1), while the first row is estimation of dx, and second row is estimation of dy. The estimation of kth factor is in the (k+1)th column of d_est

function [u_est, V_est, W_est, d_est] = Jisst_multi(X, Y, u0, lambda, tol, max_iter, rank_max, def, varargin)
    
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
        rank{1} = rx(k); rank{2} = ry(k);
        r_max = 0;
        [hat_u, hat_V, hat_W, d_x, d_y, ~, ~] = Jisst_single(X_est{k}, Y_est{k}, u_est{k}, lambda(k), tol, max_iter, r_max, rank);
        
        % update of tensor factors
        u_est{k+1} = hat_u;
        V_est{k+1} = hat_V;
        W_est{k+1} = hat_W;
        d_est(1, k+1) = d_x; % estimation of d_{x_{k}}
        d_est(2, k+1) = d_y; % estimation of d_{y_{k}}

        % deflate
        if def == 0 % subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - d_y*squeeze(ttt(tensor(hat_W*hat_W'), tensor(hat_u)));
        elseif def == 1 % project deflation
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = ttm(Y_est{k}, {Wk, Wk, uk}, [1, 2, 3]);
        elseif def == 2 % orthogonal joint factor after subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - d_y*squeeze(ttt(tensor(hat_W*hat_W'), tensor(hat_u)));
            
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');

            X_est{k+1} = ttm(X_est{k+1}, {uk}, 3);
            Y_est{k+1} = ttm(Y_est{k+1}, {uk}, 3);
        elseif def == 3 % orthogonal individual factors after subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - d_y*squeeze(ttt(tensor(hat_W*hat_W'), tensor(hat_u)));
            
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            
            X_est{k+1} = ttm(X_est{k+1}, {Vk, Vk}, [1, 2]);
            Y_est{k+1} = ttm(Y_est{k+1}, {Wk, Wk}, [1, 2]);
        end
        
        % finish the while loop
        k = k+1;
    end

    % by this function, the estimation of kth factor is (k+1)th element in
    % each cell, i.e. \hat{u}_{k} = u_est{k+1} and also for V, W
end

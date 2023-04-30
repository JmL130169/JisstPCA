% function of multi-factor matrix-tensor factorization for semi-symmetric
% tensor and rank K matrix. Using subtract deflation with imt_single

% Note: same as Jisst_multi and dJisst_multi, the estimation of kth factor
% by imt_multi is (k+1)th element of the output cells

% input of imt_multi:
% X: multi-factor semi-symmetric tensor with number of layers equals to K
% Y: rank-K matrix
% all the other agruements are the same as that in Jisst_multi

% output of imt_multi:
% estimation of each factor

function [u_est, V_est, w_est, d_est] = imt_multi(X, Y, u0, rx, lambda, tol, max_iter, def)
    
    sz = size(rx, 2); % number of layers for multi-factor iSST-PCA model
    X_est = cell(sz + 1, 1);
    Y_est = cell(sz + 1, 1);
    u_est = cell(sz + 1, 1);
    V_est = cell(sz + 1, 1);
    w_est = cell(sz + 1, 1);
    d_est = zeros(2, sz+1);

    u_est{1} = u0;
    X_est{1} = X;
    Y_est{1} = Y;

    k = 1;
    while k < sz + 1
        [hat_u, hat_V, hat_w, ~, ~, d_x, d_y] = imt_single(X_est{k}, Y_est{k}, u_est{k}, rx(k), lambda(k), tol, max_iter);
        
        % update of tensor factors
        u_est{k+1} = hat_u;
        V_est{k+1} = hat_V;
        w_est{k+1} = hat_w;
        d_est(1, k+1) = d_x;
        d_est(2, k+1) = d_y;

        % deflation
        if def == 0 % subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - d_y*hat_w*hat_u';
        else if def == 1 % project deflation
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            wk = double(eye(sz_Y(1)) - w_est{k+1}*w_est{k+1}');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = wk*Y_est{k}*uk;
        else if def == 2 % subtract deflation for X and project deflation for Y
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');            
            wk = double(eye(sz_Y(1)) - w_est{k+1}*w_est{k+1}');            
            Y_est{k+1} = wk*Y_est{k}*uk;
        else if def == 3 % project deflation for X and subtract deflation for Y
            sz_X = size(X);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = Y_est{k} - d_y*hat_w*hat_u';
        end
        
        % finish the while loop
        k = k+1;
    end
end

% Function of single-factor JisstPCA with diagonal matrix D as signal; 
% Difference from dJisst_single: estimates in each iteration are returned 

% Input of dJisst_single:

% X, Y: two single-factor semi-symmetric tensors of dimension p-p-N and q-q-N
% u0: initialization. 
% rx, ry: rank of X and Y
% lambda: scaler that represents the relative importance of each tensor
%  max_iter: maximum iteration number


% Output of dJisst_single:

% hat_u, hat_V, hat_W: estimation of factors u, V, W
% Dx, Dy: estimation of Dx and Dy
% loss: weighted Frobenious norm approximation errors

function [hat_u, hat_V, hat_W, hat_Dx, hat_Dy, loss] = dJisst_single_iter(X, Y, u0, rx, ry, lambda, max_iter, returnloss)
    
    hat_u = cell(max_iter + 1, 1); hat_V = cell(max_iter, 1); hat_W = cell(max_iter, 1);
    hat_Dx = cell(max_iter, 1); hat_Dy = cell(max_iter, 1);
    if returnloss
        loss = zeros(max_iter, 1);
    else
        loss = NaN;
    end
    hat_u{1} = u0; % initialization of u
    k = 1;

    while k < (max_iter+1)

        % update of V
        X_mat = double(ttv(X, hat_u{k}, 3));
        [hat_V{k}, ~, ~] = svds(X_mat, rx);
    
        % update of W
        Y_mat = double(ttv(Y, hat_u{k}, 3));
        [hat_W{k}, ~, ~] = svds(Y_mat, ry);
        
        % update of Dx, unlike original algorithm, we need to update Dx
        % (and Dy) at each iteration in the modified model and algorithm
        hat_Dx{k} = diag(diag(hat_V{k}'*X_mat*hat_V{k}));

        % update of Dy, similar as the update of Dx
        hat_Dy{k} = diag(diag(hat_W{k}'*Y_mat*hat_W{k}));
    
        % update of u, gtr_prod is another function which is trace product
        % of X(Y), hat_V(hat_W) and Dx(Dy). The function of gtr_prod is in
        % another .m file
        u_unnorm = lambda*gtr_prod(X, hat_V{k}, hat_Dx{k}) + (1-lambda)*gtr_prod(Y, hat_W{k}, hat_Dy{k});
        hat_u{k + 1} = u_unnorm/norm(u_unnorm);

        % estimate X and Y
        if returnloss
            hat_X = squeeze(ttt(tensor(hat_V{k}*hat_Dx{k}*hat_V{k}'), tensor(hat_u{k + 1})));
            hat_Y = squeeze(ttt(tensor(hat_W{k}*hat_Dy{k}*hat_W{k}'), tensor(hat_u{k + 1})));
            loss(k) = norm(hat_X - X)^2 * lambda + norm(hat_Y - Y)^2 * (1 - lambda);
        end
        % next iteration
        k = k+1;
    end
end





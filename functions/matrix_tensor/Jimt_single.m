% Function of single-factor matrix_tensor factorization for rank one 
% matrix and smei-symmetric tensor

% input of imt_single:

% X: single factor semi-symmetric tensor in the form of dx*VV' \circ u
% Y: rank one matrix in the form of dy*wu'
% r: rank of X
% lambda: scaler that represents the relative importance of each tensor
% tol, max_iter: tolerance value and maximum iteration number

% output of imt_single:

% estimation of each factor

function [hat_u, hat_V, hat_w, hat_X, hat_Y, d_x, d_y] = Jimt_single(X, Y, u0, r, lambda, tol, max_iter)

    % initialization
    hat_u = u0;
    k = 0;

    % update of each iteration
    while k < max_iter

        % update of V
        [a1, ~, ~] = svds(double(ttv(X, hat_u, 3)), r);
        hat_V = a1;

        % update of w
        hat_w = (Y*hat_u)/norm(Y*hat_u);

        % update of d_x and d_y
        d_x = ttt(X, squeeze(ttt(tensor(hat_V * hat_V'), tensor(hat_u))), [1, 2, 3])/r;
        d_y = hat_w'*Y*hat_u;

        % update of u
        hat_u = (lambda*tr_prod(X, hat_V)/d_x + (1-lambda)*Y'*hat_w/d_y)/norm(lambda*tr_prod(X, hat_V)/d_x + (1-lambda)*Y'*hat_w/d_y);

        % reconstruction of X and Y
        hat_X = d_x*squeeze(ttt(tensor(hat_V * hat_V'), tensor(hat_u)));
        hat_Y = d_y*hat_w*hat_u';

        k = k+1;

        % stopping criteria
        if norm(hat_X-X)/norm(X)+norm(hat_Y-Y)/norm(Y)<tol
            break
        end
    end
end

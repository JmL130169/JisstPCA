% Function of single-factor JisstPCA with diagonal matrix D as signal 

% Input of dJisst_single:

% X, Y: two single-factor semi-symmetric tensors of dimension p-p-N and q-q-N
% u0: initialization. 
% rx, ry: rank of X and Y
% lambda: scaler that represents the relative importance of each tensor
% tol, max_iter: tolerance value and maximum iteration number


% Output of dJisst_single:

% hat_u, hat_V, hat_W: estimation of factors u, V, W
% Dx, Dy: estimation of Dx and Dy
% hat_X, hat_Y: reconstruction of true parameter tensors X^{*} and Y^{*}

function [hat_u, hat_V, hat_W, Dx, Dy, hat_X, hat_Y] = dJisst_single(X, Y, u0, rx, ry, lambda, tol, max_iter)
    
    hat_u = u0; % initialization of u
    k = 1;

    while k < (max_iter+1)

        % update of V, the leading r_x singular vectors of third-mode
        % tensor multiplication with hat_u at each iteration
        X_mat = double(ttv(X, hat_u, 3));
        [ax, ~, ~] = svds(X_mat, rx);
        hat_V = ax;

        % update of W, the leading r_y singular vectors of third-mode
        % multiplication with hat_u at each iteration
        Y_mat = double(ttv(Y, hat_u, 3));
        [ay, ~, ~] = svds(Y_mat, ry);
        hat_W = ay;

        % update of Dx, unlike original algorithm, we need to update Dx
        % (and Dy) at each iteration in the modified model and algorithm
        Dx = diag(diag(hat_V'*X_mat*hat_V));

        % update of Dy, similar as the update of Dx
        Dy = diag(diag(hat_W'*Y_mat*hat_W));

        % update of u, gtr_prod is another function which is trace product
        % of X(Y), hat_V(hat_W) and Dx(Dy). The function of gtr_prod is in
        % another .m file
        u_unnorm = lambda*gtr_prod(X, hat_V, Dx) + (1-lambda)*gtr_prod(Y, hat_W, Dy);
        hat_u = u_unnorm/norm(u_unnorm);

        % estimate X and Y
        if k > 1
            hat_X_old = hat_X; hat_Y_old = hat_Y;
        end
        hat_X = squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
        hat_Y = squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));

        % next iteration
        k = k+1;

        % stopping criteria
        if k > 2 && norm(hat_X-hat_X_old)/norm(hat_X_old)+norm(hat_Y-hat_Y_old)/norm(hat_Y_old)<tol
            break
        end
    end
end








% function of JISST-PCA with diagonal matrix D as signal

% Input of dJisst_single is:
% X, Y: two single-factor semi-symmetric tensors of dimension p-p-N and q-q-N
% u0: initialization
% rx, ry: rank of X and Y
% lambda: scaler that represents the relative importance oif each tensor
% tol, max_iter: tolerance value and maximum iteration number

% Output of dJisst_single is:
% hat_u, hat_V, hat_W: estimation of factors u, V, W
% Dx, Dy: estimation of Dx and Dy
% hat_X, hat_Y: reconstruction of true parameter tensors X^{*} and Y^{*}

function [hat_u, hat_V, hat_W, Dx, Dy, hat_X, hat_Y] = dJisst_single(X, Y, u0, rx, ry, lambda, tol, max_iter)

    hat_u = u0; % initialization of u
    k = 0;

    while k < max_iter

        % update of V, the leading r_x singular vectors of third-mode
        % tensor multiplication with hat_u at each iteration
        [ax, bx, ~] = svd(double(ttv(X, hat_u, 3)));
        [~, ind_x] = sort(diag(bx));
        ax = ax(:, ind_x);
        ax = fliplr(ax);
        hat_V = ax(:, 1:rx);

        % update of W, the leading r_y singular vectors of third-mode
        % multiplication with hat_u at each iteration
        [ay, by, ~] = svd(double(ttv(Y, hat_u, 3)));
        [~, ind_y] = sort(diag(by));
        ay = ay(:, ind_y);
        ay = fliplr(ay);
        hat_W = ay(:, 1:ry);

        % update of Dx, unlike original algorithm, we need to update Dx
        % (and Dy) at each iteration in the modified model and algorithm
        Dx = diag(diag(hat_V'*(double(ttv(X, hat_u, 3)))*hat_V));

        % update of Dy, similar as the update of Dx
        Dy = diag(diag(hat_W'*(double(ttv(Y, hat_u, 3)))*hat_W));

        % update of u, gtr_prod is another function which is trace product
        % of X(Y), hat_V(hat_W) and Dx(Dy). The function of gtr_prod is in
        % another .m file
        u_unnorm = lambda*gtr_prod(X, hat_V, Dx) + (1-lambda)*gtr_prod(Y, hat_W, Dy);
        hat_u = u_unnorm/norm(u_unnorm);

        % estimate X and Y
        hat_X = squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
        hat_Y = squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));

        % next iteration
        k = k+1;

        % stopping criteria
        if norm(hat_X-X)/norm(X)+norm(hat_Y-Y)/norm(Y)<tol
            break
        end
    end
end








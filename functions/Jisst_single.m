% This function is for single-factor JisstPCA Algorithm 

% Input of Jisst_single:

% X, Y: two single-factor semi-symmetric tensors of dimension p-p-N and q-q-N
% u0: initialization. 
% rx, ry: rank of X and Y
% lambda: scaler that represents the relative importance of each tensor
% tol, max_iter: tolerance value and maximum iteration number


% Output of Jisst_single:

% hat_u, hat_V, hat_W: estimation of factors u, V, W
% dx, dy: estimation of dx and dy
% hat_X, hat_Y: reconstruction of true parameter tensors X^{*} and Y^{*}, 
% by hat_X = dx*hat_V hat_V' \circ hat_u and so is for hat_Y

function [hat_u, hat_V, hat_W, dx, dy, hat_X, hat_Y] = Jisst_single(X, Y, u0, rx, ry, lambda, tol, max_iter)

    % Note: tensor_toolbox by Kolda is needed for tensor algebra in this
    % function
    %
    % Jisst_single computes two integrated semi-symmetric tensors decomposition,
    % where the input X, Y are p-by-p-by-N and q-by-q-by-N semi-symmetric
    % tensors, which can be decomposed in terms of "True semi-symmetric
    % tensor + semi-symmetric noise".
    %
    % Suppose "X = dx V V'*u + noise_X", "Y = dy W W'*u + noise_Y". V, W are
    % p-by-rx, q-by-ry orthogonal matrix, u is a unit vector of dimension N
    % which is the joint structure. Here * stands for outer product of
    % matrix and vector.
    %
    % The input X, Y are the observed semi-symmetric tensor, u0 is
    % initialization of joint tensor factor u, rx and ry are
    % rank of V and W. Lambda is a data-driven scaler to
    % describe the relative importance of the observations.
    %
    % Note: The data type of input X, Y should be 'tensor'.


    k = 1;
    hat_u = u0; % at the 0th iteration, estimation of joint tensor factor u
    % is just the initialization value of u;
    
    % update of each tensor factor
    while k < (max_iter + 1)
    
        % update of V
        % at each iteration, the update of V is the rx leading singular vectors
        % of X \times_3 u^{-}, the tensor product of X and the estimation
        % of joint tensor factor u from the previous step.
        %
        % Note: at each iterstion, the estimate of V should be an orthogonal
        % matrix. Also, svd function in MATLAB have already arrange the
        % singular value/vector in an decreasing order, we arrange them again
        % to make sure we get the leading singular vectors
        X_mat = double(ttv(X, hat_u, 3));
        [hat_V, ~, ~] = svds(X_mat, rx);
    
        % update of W
        % at each iteration, the update of W is the ry leading eigenvectors
        % of Y \times_3 u^{-}, the tensor product of Y and the estimation
        % of joint tensor factor u from the previous step.
        %
        % Note: at each iterstion, the estimate of W should be an orthogonal
        % matrix.
        Y_mat = double(ttv(Y, hat_u, 3));
        [hat_W, ~, ~] = svds(Y_mat, ry);
    
        % update of u
        % the update of u is a weighted summation of two trace products:
        % trace product of X and V (V is from the current step), trace
        % product of Y and W (W is from the current step). And normalize
        % this vector given that u should be a unit vector.
        %
        % Note: the function of trace product is tr_prod given by another
        % file in this toolbox.
        u_unnorm = lambda*tr_prod(X, hat_V)+(1-lambda)*tr_prod(Y, hat_W);
        hat_u = u_unnorm/norm(u_unnorm);

        % estimate X and Y
        % after obtaining the estimation of u, V, W, since dx*rx equals to
        % the inner product of X and itself (so is dy*ry), the estimation of
        % tensor signal dx is <X, hat_V hat_V'*hat_u>/rx (so is dy).
        %
        % then the estimation of the true X can be obtained in terms of
        % hat_V, hat_u, dx (so is the estimation of the true Y).
        if k > 1
            hat_X_old = hat_X; hat_Y_old = hat_Y;
        end
        dx = trace(hat_V'*X_mat*hat_V) / rx;% the estimated signal of X;
        dy = trace(hat_W'*Y_mat*hat_W) / ry;% the estimated signal of Y;

        hat_X = dx * squeeze(ttt(tensor(hat_V * hat_V'), tensor(hat_u))); % estimation of X;
        hat_Y = dy * squeeze(ttt(tensor(hat_W * hat_W'), tensor(hat_u))); % estimation of Y;

        % next iteration
        k = k+1;

        % stopping criteria
        if k > 2 && norm(hat_X-hat_X_old)/norm(hat_X_old)+norm(hat_Y-hat_Y_old)/norm(hat_Y_old)<tol
            break
        end
    end
end

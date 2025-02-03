% This function is for single-factor sstPCA Algorithm 

%%%%%%%%%%%%%%%%%%%%%%%%
% Input of sst_single: %
%%%%%%%%%%%%%%%%%%%%%%%%

% X: Single-factor semi-symmetric tensors of dimension p-p-N
% u0: initialization. 
% rx: rank of X
% tol, max_iter: tolerance value and maximum iteration number

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output of Jisst_single: %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% hat_u, hat_V: estimation of factors u, V
% dx: estimation of dx
% hat_X: reconstruction of true parameter tensors X^{*}, 
% by hat_X = dx*hat_V hat_V' \circ hat_u and so is for hat_Y

function [hat_u, hat_V, dx, hat_X] = sst_single(X, u0, rx, tol, max_iter)

    % Note: tensor_toolbox by Kolda is needed for tensor algebra in this
    % function
    %
    % sst_single computes semi-symmetric tensor decomposition in Weylandt & Michailidis (2022),
    % where the input X is a p-by-p-by-N semi-symmetric tensor, 
    % which can be decomposed in terms of "True semi-symmetric tensor + semi-symmetric noise".
    %
    % Suppose "X = dx V V'*u + noise_X". V is a
    % p-by-rx orthogonal matrix, u is a unit vector of dimension N
    % which is the joint structure. Here * stands for outer product of
    % matrix and vector.
    %
    % The input X is the observed semi-symmetric tensor, u0 is
    % initialization of joint tensor factor u, rx is rank of V.
    %
    % Note: The data type of input X should be 'tensor'.


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
    
        % update of u
        % the update of u is a weighted summation of two trace products:
        % trace product of X and V (V is from the current step), trace
        % product of Y and W (W is from the current step). And normalize
        % this vector given that u should be a unit vector.
        %
        % Note: the function of trace product is tr_prod given by another
        % file in this toolbox.
        u_unnorm = tr_prod(X, hat_V);
        hat_u = u_unnorm/norm(u_unnorm);

        % estimate X
        % after obtaining the estimation of u, V, since dx*rx equals to
        % the inner product of X and itself, the estimation of
        % tensor signal dx is <X, hat_V hat_V'*hat_u>/rx.
        %
        % then the estimation of the true X can be obtained in terms of
        % hat_V, hat_u, dx.
        if k > 1
            hat_X_old = hat_X;
        end
        dx = trace(hat_V'*X_mat*hat_V) / rx;% the estimated signal of X;

        hat_X = dx * squeeze(ttt(tensor(hat_V * hat_V'), tensor(hat_u))); % estimation of X;

        % next iteration
        k = k+1;

        % stopping criteria
        if k > 2 && norm(hat_X-hat_X_old)/norm(hat_X_old)<tol
            break
        end
    end
end
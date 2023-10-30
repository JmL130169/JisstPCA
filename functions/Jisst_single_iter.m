% This function is for single-factor JisstPCA Algorithm 
% Difference from Jisst_single: estimates in each iteration are returned 

% Input of Jisst_single:

% X, Y: two single-factor semi-symmetric tensors of dimension p-p-N and q-q-N
% u0: initialization. 
% rx, ry: rank of X and Y
% lambda: scaler that represents the relative importance of each tensor
% max_iter: maximum iteration number


% Output of Jisst_single:

% hat_u, hat_V, hat_W: estimation of factors u, V, W in each iteration
% loss: weighted Frobenious norm approximation errors

function [hat_u, hat_V, hat_W, loss] = Jisst_single_iter(X, Y, u0, rx, ry, lambda, max_iter, returnloss)

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
    
    hat_u = cell(max_iter + 1, 1); hat_V = cell(max_iter, 1); hat_W = cell(max_iter, 1);
    if returnloss
        loss = zeros(max_iter, 1);
    else
        loss = NaN;
    end
    k = 1;
    hat_u{1} = u0; % at the 0th iteration, estimation of joint tensor factor u
    % is just the initialization value of u;
    % update of each tensor factor
    while k < (max_iter + 1)
    
        % update of V
        [hat_V{k}, ~, ~] = svds(double(ttv(X, hat_u{k}, 3)), rx);
    
        % update of W
        [hat_W{k}, ~, ~] = svds(double(ttv(Y, hat_u{k}, 3)), ry);
    
        % update of u
        hat_u{k + 1} = (lambda*tr_prod(X, hat_V{k})+(1-lambda)*tr_prod(Y, hat_W{k}))/norm(lambda*tr_prod(X, hat_V{k})+(1-lambda)*tr_prod(Y, hat_W{k}));

        % estimate X and Y
        if returnloss
            hat_X0 = squeeze(ttt(tensor(hat_V{k} * hat_V{k}'), tensor(hat_u{k + 1})));
            hat_Y0 = squeeze(ttt(tensor(hat_W{k} * hat_W{k}'), tensor(hat_u{k + 1})));
            dx = ttt(X, hat_X0, 1:3)/rx; % the estimated signal of X;
            dy = ttt(Y, hat_Y0, 1:3)/ry; % the estimated signal of Y;
            hat_X = tensor(dx*hat_X0); % estimation of X;
            hat_Y = tensor(dy*hat_Y0); % estimation of Y;
            loss(k) = norm(hat_X - X)^2 * lambda + norm(hat_Y - Y)^2 * (1 - lambda);
        end
        % next iteration
        k = k + 1;
    end

end
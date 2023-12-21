% This function uses alternating minimization to solve single-factor JisstPCA
% Input:
% X, Y: two single-factor semi-symmetric tensors of dimension p-p-N and q-q-N
% u0: initialization. 
% rx, ry: rank of X and Y
% lambda: scaler that represents the relative importance of each tensor
% max_iter: maximum iteration number


% Output:

% hat_u, hat_V, hat_W, hat_d_x, hat_d_y: estimation of factors u, V, W, d_x, d_y in each iteration
% loss: weighted Frobenious norm approximation errors

function [hat_u, hat_V, hat_W, hat_d_x, hat_d_y, loss] = Jisst_single_alm(X, Y, u0, rx, ry, lambda, max_iter, returnloss)

    % Note: tensor_toolbox by Kolda is needed for tensor algebra in this
    % function
    %
    % Jisst_single_alm computes two integrated semi-symmetric tensors decomposition,
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
    hat_d_x = zeros(max_iter, 1); hat_d_y = zeros(max_iter, 1);
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
        Xu =  double(ttv(X, hat_u{k}, 3));
        [hat_V{k}, ~, ~] = svds(Xu, rx);
    
        % update of W
        Yu =  double(ttv(Y, hat_u{k}, 3));
        [hat_W{k}, ~, ~] = svds(Yu, ry);
        
        % update of d_x
        hat_d_x(k) = trace(hat_V{k}' * Xu * hat_V{k})/rx;

        % update of d_y
        hat_d_y(k) = trace(hat_W{k}' * Yu * hat_W{k})/ry;

        % update of u
        hat_u{k + 1} = (lambda*hat_d_x(k)*tr_prod(X, hat_V{k})+(1-lambda)*hat_d_y(k)*tr_prod(Y, hat_W{k}));
        hat_u{k + 1} = hat_u{k + 1}/norm(hat_u{k + 1});

        % estimate X and Y
        if returnloss
            hat_X0 = squeeze(ttt(tensor(hat_V{k} * hat_V{k}'), tensor(hat_u{k + 1})));
            hat_Y0 = squeeze(ttt(tensor(hat_W{k} * hat_W{k}'), tensor(hat_u{k + 1})));
            hat_X = tensor(hat_d_x(k)*hat_X0); % estimation of X;
            hat_Y = tensor(hat_d_y(k)*hat_Y0); % estimation of Y;
            loss(k) = norm(hat_X - X)^2 * lambda + norm(hat_Y - Y)^2 * (1 - lambda);
        end
        % next iteration
        k = k + 1;
    end

end
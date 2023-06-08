% Function of the joint-integrated matrix_tensor PCA 

% Input of JimtPCA:

% (a). Required input parameters:

% X: semi-symmetric tensor of dimension p-p-N, with K layers
% Y: matrix of dimension q-N, of rank K
% K: number of layers, defined by the user

% (b). Optional input parameters (if not identified, JisstPCA will automatically
% apply default ones):

% u0: initialization. The default value is spectral initialization
% rx: rank of X, which is vector of dimension K. The default
% value is obtained by BIC + deflation scheme
% lambda: scaler for each layer, which is a vector of dimension K. The
% default value is 0.5 for each element
% tol, max_iter: tolerance value and maximum iteration number. The default
% values are tol = 0.0001 and max_iter = 20 
% deflation: deflation = 0 is subtract deflation, deflation = 1 is project 
% deflation, deflation = 2 is project deflation for only u after subtract 
% deflation, def = 3 is project deflation for only V and W after subtract
% deflation. The default deflation strategy is deflation = 0, which is
% subtract deflation

% (c). Input parameters used in BIC for selecting the best rx:
% rank_max: if the true rank rx is unknown, rank_max is the largest
% possible that will be tried using BIC method. The default value is
% rank_max = 5


% Output of JisstPCA:

% u_est, V_est, W_est: cells with K elements. The estimation of kth factor 
% of u, V, W are the kth factor in u_est, V_est, W_est. For example, the
% estimation of uk is u_est{k}
% d_est: matrix of dimension 2*K, while the first row is estimation of dx, 
% and second row is estimation of dy. For example, the estimation of dx_k
% is d_est(1, k)

function [u_est, V_est, w_est, d_est] = JimtPCA(X, Y, K, varargin)

    % create an input parser object
    p = inputParser;

    % add required and optional parameters
    addRequired(p, 'X');
    addRequired(p, 'Y');
    addRequired(p, 'K');
    addParameter(p, 'u0', NaN);
    addParameter(p, 'rx', NaN);
    addParameter(p, 'lambda', NaN);
    addParameter(p, 'tol', NaN);
    addParameter(p, 'max_iter', NaN);
    addParameter(p, 'deflation', NaN);
    addParameter(p, 'rank_max', NaN);

    % parse the inputs
    parse(p, X, Y, K, varargin{:});

    % get the values
    X = p.Results.X;
    Y = p.Results.Y;
    K = p.Results.K;
    u0 = p.Results.u0;
    rx = p.Results.rx;
    lambda = p.Results.lambda;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
    deflation = p.Results.deflation;
    rank_max = p.Results.rank_max;
    
    % check if the user proveides lambda, tol, max_iter and
    % deflation
    if isnan(lambda)
        lam = 0.5;
        lambda = lam*ones(K, 1);
    end
    if isnan(tol)
        tol = 0.0001;
    end
    if isnan(max_iter)
        max_iter = 20;  
    end
    if isnan(deflation)
        deflation = 0; 
    end
    if isnan(rank_max)
        rank_max = 5;  
    end

    % default value for u0, rx and ry
    if isnan(u0)
        u0 = init_mt(X, Y); 
    end
    if isnan(rx)
        rx = bic_mt(X, Y, rank_max, K, u0, lambda, tol, max_iter, deflation);
    end
    
    %% main function when all the hyperparameters are known

    sz = size(rx, 2); 
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
    while k < K + 1
        [hat_u, hat_V, hat_w, ~, ~, d_x, d_y] = Jimt_single(X_est{k}, Y_est{k}, u_est{k}, rx(k), lambda(k), tol, max_iter);
        
        % update of tensor factors
        u_est{k+1} = hat_u;
        V_est{k+1} = hat_V;
        w_est{k+1} = hat_w;
        d_est(1, k+1) = d_x;
        d_est(2, k+1) = d_y;

        % different deflation strategy

        if deflation == 0 % subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - d_y*hat_w*hat_u';
        elseif deflation == 1 % project deflation
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            wk = double(eye(sz_Y(1)) - hat_w*hat_w');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = wk*Y_est{k}*uk;
        elseif deflation == 2 % project defaltion on joint factor
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - d_y*hat_w*hat_u';
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            wk = double(eye(sz_Y(1)) - hat_w*hat_w');
            X_est{k+1} = ttm(X_est{k+1}, {uk}, 3);
            Y_est{k+1} = Y_est{k+1}*uk;
        elseif deflation == 3 % project deflation on individual factor
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - d_y*hat_w*hat_u';
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - hat_u*hat_u');
            Vk = double(eye(sz_X(1)) - hat_V*hat_V');
            wk = double(eye(sz_Y(1)) - hat_w*hat_w');
            X_est{k+1} = ttm(X_est{k+1}, {Vk, Vk}, [1, 2]);
            Y_est{k+1} = wk*Y_est{k+1};
        end

        % finish the while loop
        k = k+1;
    end

    % make index consistent
    u_est(1) = [];
    V_est(1) = [];
    w_est(1) = [];
    d_est(:, 1) = [];
end



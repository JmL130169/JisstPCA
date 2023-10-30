% Function of the model joint-integrated semi-symmetric tensor PCA, 
% diagonal signal D

% Input of dJisstPCA:

% (a). Required input parameters:

% X, Y: two multi-factor semi-symmetric tensors of dimension p-p-N and q-q-N, with K layers
% K: number of layers, defined by the user

% (b). Optional input parameters (if not identified, JisstPCA will automatically
% apply default ones):

% u0: initialization for all layers. If not given by the user, spectral initialization will be used for each
% layer. 
% rx, ry: rank of X and Y, both are vectors of dimension K. The default
% value is obtained by BIC + deflation scheme
% lambda: scaler for each layer, which is a vector of dimension K. The
% default value is norm(X)/(norm(X)+norm(Y)) for each element
% tol, max_iter: tolerance value and maximum iteration number. The default
% values are tol = 0.00001 and max_iter = 20 
% deflation: deflation = 0 is subtract deflation, deflation = 1 is project 
% deflation, deflation = 2 is project deflation for only u after subtract 
% deflation, def = 3 is project deflation for only V and W after subtract
% deflation. The default deflation strategy is deflation = 0, which is
% subtract deflation

% (c). Input parameters used in BIC for selecting the best rx, ry:
% rank_max: if the true rank rx, ry are unknown, rank_max is the largest
% possible that will be tried using BIC method. The default value is
% rank_max = 5
% method: when using BIC for selecting ranks, method = 1 means the user
% does not assume X and Y are of the same rank, while method = 2 means the
% user would assume X and Y are of the same rank. The default value is
% method = 1


% Output of dJisstPCA:

% u_est, V_est, W_est: cells with K elements. The estimation of kth factor 
% of u, V, W are the kth factor in u_est, V_est, W_est. For example, the
% estimation of uk is u_est{k}
% Dx_est, Dy_est; cells with K elements. The estimation of kth signal of X
% (Dx_k) or Y (Dy_k) are the kth elements in Dx_est or Dy_est

function [u_est, V_est, W_est, Dx_est, Dy_est] = dJisstPCA(X, Y, K, varargin)

    % create an input parser object
    p = inputParser;

    % add required and optional parameters
    addRequired(p, 'X');
    addRequired(p, 'Y');
    addRequired(p, 'K');
    addParameter(p, 'u0', NaN);
    addParameter(p, 'rx', NaN);
    addParameter(p, 'ry', NaN);
    addParameter(p, 'lambda', NaN);
    addParameter(p, 'tol', NaN);
    addParameter(p, 'max_iter', NaN);
    addParameter(p, 'deflation', NaN);
    addParameter(p, 'rank_max', NaN);
    addParameter(p, 'method', NaN);

    % parse the inputs
    parse(p, X, Y, K, varargin{:});

    % get the values
    X = p.Results.X;
    Y = p.Results.Y;
    K = p.Results.K;
    u0 = p.Results.u0;
    rx = p.Results.rx;
    ry = p.Results.ry;
    lambda = p.Results.lambda;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
    deflation = p.Results.deflation;
    rank_max = p.Results.rank_max;
    method = p.Results.method;

    % check if the user proveides lambda, tol, max_iter and
    % deflation
    if isnan(lambda)
        lam = norm(X)/(norm(X)+norm(Y));
        lambda = lam*ones(K, 1);
    end
    if isnan(tol)
        tol = 0.00001;
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
    if isnan(method)
        method = 1; 
    end

    % default value for u0, rx and ry
    if isnan(u0)
        spectral_init = true;
        u0 = init(X, Y, lambda(1)); 
    else
        spectral_init = false;
    end
    if max(isnan(rx))||max(isnan(ry))
        [rx_tmp, ry_tmp] = bic_diag(X, Y, rank_max, K, u0, spectral_init, lambda, tol, max_iter, method, deflation);
        if max(isnan(rx))
            rx = rx_tmp;
        end
        if max(isnan(ry))
            ry = ry_tmp;
        end
    end
    
    %% main function when all the hyperparameters are known
    sz = size(rx, 2);
    X_est = cell(sz + 1, 1);
    Y_est = cell(sz + 1, 1);
    u_est = cell(sz + 1, 1);
    V_est = cell(sz + 1, 1);
    W_est = cell(sz + 1, 1);
    Dx_est = cell(sz + 1, 1);
    Dy_est = cell(sz + 1, 1);

    % set initialization and let X^{1} = X, Y^{1} = Y as in our algorithm
    u_est{1} = u0;
    X_est{1} = X;
    Y_est{1} = Y;

    k = 1;
    while k < (K + 1)
        if (k > 1) && spectral_init
            u0 = init(X_est{k}, Y_est{k}, lambda(k)); 
        end
        [hat_u, hat_V, hat_W, Dx, Dy, ~, ~] = dJisst_single(X_est{k}, Y_est{k}, u0, rx(k), ry(k), lambda(k), tol, max_iter);
        
        % update of tensor factors
        u_est{k+1} = hat_u;
        V_est{k+1} = hat_V;
        W_est{k+1} = hat_W;
        Dx_est{k+1} = Dx;
        Dy_est{k+1} = Dy;

        % deflation
        if deflation == 0 % subtract deflation
            X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));
        elseif deflation == 1 % project deflation
            sz_X = size(X);
            sz_Y = size(Y);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
            Y_est{k+1} = ttm(Y_est{k}, {Wk, Wk, uk}, [1, 2, 3]);
        elseif deflation == 2 % orthogonal joint factor after subtract deflation
            X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));
            
            sz_X = size(X);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            
            X_est{k+1} = ttm(X_est{k+1}, {uk}, 3);
            Y_est{k+1} = ttm(Y_est{k+1}, {uk}, 3);
        elseif deflation == 3 % orthogonal individual factors after subtract deflation
            X_est{k+1} = X_est{k} - squeeze(ttt(tensor(hat_V*Dx*hat_V'), tensor(hat_u)));
            Y_est{k+1} = Y_est{k} - squeeze(ttt(tensor(hat_W*Dy*hat_W'), tensor(hat_u)));
            
            sz_X = size(X);
            sz_Y = size(Y);
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            Wk = double(eye(sz_Y(1)) - W_est{k+1}*W_est{k+1}');
            
            X_est{k+1} = ttm(X_est{k+1}, {Vk, Vk}, [1, 2]);
            Y_est{k+1} = ttm(Y_est{k+1}, {Wk, Wk}, [1, 2]);
        end
        
        % next iteration
        k = k+1;
    end

    % make index consistent
    u_est(1) = [];
    V_est(1) = [];
    W_est(1) = [];
    Dx_est(1) = [];
    Dy_est(1) = [];

end








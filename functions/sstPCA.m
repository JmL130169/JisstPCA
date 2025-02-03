% Main function of the model semi-symmetric tensor PCA W & M (2022)

%%%%%%%%%%%%%%%%%%%%
% Input of sstPCA: %
%%%%%%%%%%%%%%%%%%%%

% (a). Required input parameters:

% X: multi-factor semi-symmetric tensor of dimension p-p-N, with K layers
% K: number of layers, defined by the user

% (b). Optional input parameters (if not identified, sstPCA will automatically
% apply default ones):

% u0: initialization. The default value is spectral initialization
% rx: rank of X, that is vector of dimension K. The default
% value is obtained by BIC + deflation scheme
% tol, max_iter: tolerance value and maximum iteration number. The default
% values are tol = 0.00001 and max_iter = 20 
% deflation: deflation = 0 is subtract deflation, deflation = 1 is project 
% deflation, deflation = 2 is project deflation for only u after subtract 
% deflation, def = 3 is project deflation for only V after subtract
% deflation. The default deflation strategy is deflation = 0, which is
% subtract deflation

% (c). Input parameters used in BIC for selecting the best rx:
% rank_max: if the true rank rx is unknown, rank_max is the largest
% possible that will be tried using BIC method. The default value is
% rank_max = 5

%%%%%%%%%%%%%%%%%%%%%%%
% Output of JisstPCA: %
%%%%%%%%%%%%%%%%%%%%%%%

% u_est, V_est: cells with K elements. The estimation of kth factor 
% of u, V are the kth factor in u_est, V_est. For example, the
% estimation of uk is u_est{k}
% d_est: matrix of dimension 1*K representing the signal of each layer

function [u_est, V_est, d_est] = sstPCA(X, K, varargin)

    % create an input parser object
    p = inputParser;

    % add required and optional parameters
    addRequired(p, 'X');
    addRequired(p, 'K');
    addParameter(p, 'u0', NaN);
    addParameter(p, 'rx', NaN);
    addParameter(p, 'tol', NaN);
    addParameter(p, 'max_iter', NaN);
    addParameter(p, 'deflation', NaN);
    addParameter(p, 'rank_max', NaN);

    % parse the inputs
    parse(p, X, K, varargin{:});

    % get the values
    X = p.Results.X;
    K = p.Results.K;
    u0 = p.Results.u0;
    rx = p.Results.rx;
    tol = p.Results.tol;
    max_iter = p.Results.max_iter;
    deflation = p.Results.deflation;
    rank_max = p.Results.rank_max;

    % check if the user provides lambda, tol, max_iter and
    % deflation
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

    % default value for u0, rx
    if isnan(u0)
        spectral_init = true;
        u0 = init_sst(X); 
    else
        spectral_init = false;
    end

    if isnan(rx)
        rx_tmp = bic_def_uni(X, rank_max, K, u0, spectral_init, tol, max_iter, deflation);
        rx = rx_tmp;
    end
    
    %% main function when all the hyperparameters are known
    sz = size(rx);
    X_est = cell(sz(2) + 1, 1);
    u_est = cell(sz(2) + 1, 1);
    V_est = cell(sz(2) + 1, 1);
    d_est = zeros(1, sz(2)+1);

    u_est{1} = u0;
    X_est{1} = X;

    k = 1;
    while k < K + 1
        if (k > 1) && spectral_init
            u0 = init_sst(X_est{k}); 
        end
        [hat_u, hat_V, d_x, ~] = sst_single(X_est{k}, u0, rx(k), tol, max_iter);
        
        % update of tensor factors
        u_est{k+1} = hat_u;
        V_est{k+1} = hat_V;
        d_est(k+1) = d_x; % estimation of d_{x_{k}}

        % deflate
        if deflation == 0 % subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
         elseif deflation == 1 % project deflation
            sz_X = size(X);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            X_est{k+1} = ttm(X_est{k}, {Vk, Vk, uk}, [1, 2, 3]);
        elseif deflation == 2 % orthogonal joint factor after subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            
            sz_X = size(X);
            uk = double(eye(sz_X(3)) - u_est{k+1}*u_est{k+1}');
            
            X_est{k+1} = ttm(X_est{k+1}, {uk}, 3);
        elseif deflation == 3 % orthogonal individual factors after subtract deflation
            X_est{k+1} = X_est{k} - d_x*squeeze(ttt(tensor(hat_V*hat_V'), tensor(hat_u)));
            
            sz_X = size(X);
            Vk = double(eye(sz_X(1)) - V_est{k+1}*V_est{k+1}');
            
            X_est{k+1} = ttm(X_est{k+1}, {Vk, Vk}, [1, 2]);
        end
        
        % finish the while loop
        k = k+1;
    end

    % make index consistent
    u_est(1) = [];
    V_est(1) = [];
    d_est(:, 1) = [];
end

%% Example
% Given X, K, if all the other hyperparameters are unknown, then we can
% just implement sstPCA(X, K)
% Given X, K, if partial (or all) hyperparameters are known, for example
% u0 and rx, then we should use this function as sstPCA(X, K, 'u0', u0, 'rx', rx) 
% to specify which arguements are available
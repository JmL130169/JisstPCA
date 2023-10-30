%% comparative study, network data
addpath('/Users/lilybean/Downloads/software/tensor_toolbox-v3.2.1') 
addpath('/Applications/Files/GitHub/JisstPCA_BIC_log/functions')
addpath('/Applications/Files/GitHub/JisstPCA_BIC_log/example')


p = 150; q_vec = [150 100]; N_vec = [200 50]; 
rx = [3 2]; ry = rx;
K = 2;

vnames = {'q', 'N', 'seed', 'method', 'factor name', 'factor number', 'err'};
results = array2table(zeros(0,7), 'VariableNames',vnames);
u = cell(K, 1); V = cell(K, 1); W = cell(K, 1);
filename = "/Applications/Files/research/JisstPCA/JisstPCA_new/network/results.csv";
q = 150; N = 200;
    % probability tensor construction
    px1 = table2array(readtable(sprintf("/Applications/Files/research/JisstPCA/JisstPCA_new/network/px1_q%d.csv", q)));
    px2 = table2array(readtable(sprintf("/Applications/Files/research/JisstPCA/JisstPCA_new/network/px2_q%d.csv", q)));
    py1 = table2array(readtable(sprintf("/Applications/Files/research/JisstPCA/JisstPCA_new/network/py1_q%d.csv", q)));
    py2 = table2array(readtable(sprintf("/Applications/Files/research/JisstPCA/JisstPCA_new/network/py2_q%d.csv", q)));
    rng(2023)
    thr = 0.75;
        ur = rand(N, 1);
        u1 = 1 * (ur < thr);
        u2 = 1 - u1;
    
        % true probability tensors
        X1 = squeeze(ttt(tensor(px1), tensor(u1)));
        X2 = squeeze(ttt(tensor(px2), tensor(u2)));
        Y1 = squeeze(ttt(tensor(py1), tensor(u1)));
        Y2 = squeeze(ttt(tensor(py2), tensor(u2)));
        X_t = X1 + X2;
        Y_t = Y1 + Y2;
    
        % ground truth factors
        px_1 = (eye(p)-ones(p, 1)*ones(p, 1)'/p)*px1*(eye(p)-ones(p, 1)*ones(p, 1)'/p);
        px_2 = (eye(p)-ones(p, 1)*ones(p, 1)'/p)*px2*(eye(p)-ones(p, 1)*ones(p, 1)'/p);
        py_1 = (eye(q)-ones(q, 1)*ones(q, 1)'/q)*py1*(eye(q)-ones(q, 1)*ones(q, 1)'/q);
        py_2 = (eye(q)-ones(q, 1)*ones(q, 1)'/q)*py2*(eye(q)-ones(q, 1)*ones(q, 1)'/q);
        [V{1},~,~] =svds(px_1, rx(1) - 1); [V{2},~,~] =svds(px_2, rx(2) - 1);
        [W{1},~,~] =svds(py_1, ry(1) - 1); [W{2},~,~] =svds(py_2, ry(2) - 1);
        u{1} = u1/norm(u1); u{2} = u2/norm(u2);
    
        % data tensors
        seed = 1;
            rng(seed)
            X = zeros(p, p, N); Y = zeros(q, q, N);X_obs = zeros(p, p, N); Y_obs = zeros(q, q, N);
            for i = 1 : N
                X_tmp = tril(double(rand(p, p) < X_t(:, :, i)), -1);
                X(:, :, i) = X_tmp + X_tmp';
                X_obs(:, :, i) = (eye(p)-ones(p, 1)*ones(p, 1)'/p)*X(:, :, i)*(eye(p)-ones(p, 1)*ones(p, 1)'/p);
                Y_tmp = tril(double(rand(q, q) < Y_t(:, :, i)), -1);
                Y(:, :, i) = Y_tmp + Y_tmp';
                Y_obs(:, :, i) = (eye(q)-ones(q, 1)*ones(q, 1)'/q)*Y(:, :, i)*(eye(q)-ones(q, 1)*ones(q, 1)'/q);      
            end
            X_obs = tensor(X_obs); Y_obs = tensor(Y_obs);
    
            vnames = {'method', 'FactorName', 'FactorNumber', 'RowInd', 'ColInd', 'Val'};
            VW_results = array2table(zeros(0, 6), 'VariableNames', vnames);
            vnames = {'method', 'FactorNumber', 'RowInd', 'Val'};
            u_results = array2table(zeros(0, 4), 'VariableNames', vnames);

            VW_results = AddStructuredValVW(VW_results, V, W, "truth");
            u_results = AddStructuredValU(u_results, u, "truth");


            deflation = 0;
            %Jisst-PCA, BIC & oracle ranks, orthogonal & subtraction deflation
            tic
            [u_est, V_est, W_est, ~] = JisstPCA(X_obs, Y_obs, 2, 'rx', rx - 1, 'ry', ry - 1,...
                'deflation', deflation, 'max_iter', 10);
            toc 
            results = add_results_structured(results, q, N, seed,...
                    "JisstPCA", find_dist(u_est, V_est, W_est, u, V, W));
            VW_results = AddStructuredValVW(VW_results, V_est, W_est, "JisstPCA");
            u_results = AddStructuredValU(u_results, u_est, "JisstPCA");        
        

            %G-JisstPCA
            tic
            [u_est, V_est, W_est, Dx_est, Dy_est] = dJisstPCA(X_obs, Y_obs, 2, 'rx', rx - 1, 'ry', ry - 1, ...
                'deflation', deflation, 'max_iter', 10);
            toc
            results = add_results_structured(results, q, N, seed,...
                    "G-JisstPCA", find_dist(u_est, V_est, W_est, u, V, W));
            VW_results = AddStructuredValVW(VW_results, V_est, W_est, 'G-JisstPCA');
            u_results = AddStructuredValU(u_results, u_est, 'G-JisstPCA'); 


            % Method 3: iHOSVD
            [u_est_mat, V_est_mat, W_est_mat, ~, ~] = iHOSVD(X_obs, Y_obs, sum(rx - 1), sum(ry - 1), 2);
            u_est = cell(K, 1); V_est = cell(K, 1); W_est = cell(K, 1); start_V = 1; start_W = 1;
            for k = 1 : K
                u_est{k} = u_est_mat(:, k);
                V_est{k} = V_est_mat(:, start_V : (start_V + rx(k) - 2));
                W_est{k} = W_est_mat(:, start_W : (start_W + ry(k) - 2));
                start_V = start_V + rx(k) - 1; start_W = start_W + ry(k) - 1;
            end
            results = add_results_structured(results, q, N, seed,...
                    "iHOSVD", find_dist(u_est, V_est, W_est, u, V, W));
            VW_results = AddStructuredValVW(VW_results, V_est, W_est, "iHOSVD");
            u_results = AddStructuredValU(u_results, u_est, "iHOSVD"); 


            % Method 4: iHOOI
            [u_est_mat, V_est_mat, W_est_mat, ~, ~] = iHOOI(X_obs, Y_obs, sum(rx - 1), sum(ry - 1), 10, 2);
            u_est = cell(K, 1); V_est = cell(K, 1); W_est = cell(K, 1); start_V = 1; start_W = 1;
            for k = 1 : K
                u_est{k} = u_est_mat(:, k);
                V_est{k} = V_est_mat(:, start_V : (start_V + rx(k) - 2));
                W_est{k} = W_est_mat(:, start_W : (start_W + ry(k) - 2));
                start_V = start_V + rx(k) - 1; start_W = start_W + ry(k) - 1;
            end
            results = add_results_structured(results, q, N, seed,...
                    "iHOOI", find_dist(u_est, V_est, W_est, u, V, W));
            VW_results = AddStructuredValVW(VW_results, V_est, W_est, "iHOOI");
            u_results = AddStructuredValU(u_results, u_est, "iHOOI");


            VW_filename = sprintf("/Applications/Files/research/JisstPCA/JisstPCA_new/network/results/factorVW_q%d_N%d_seed%d.csv", q, N, seed);
            u_filename = sprintf("/Applications/Files/research/JisstPCA/JisstPCA_new/network/results/factorU_q%d_N%d_seed%d.csv", q, N, seed);
            writetable(VW_results, VW_filename);
            writetable(u_results, u_filename);
        
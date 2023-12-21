%% comparative study, network data
% Software dependency: tensor toolbox (https://www.tensortoolbox.org)
% RI & ARI computation functions (https://www.mathworks.com/matlabcentral/fileexchange/130779-rand-and-adjusted-rand-index-calculator-for-cluster-analysis)
addpath('JisstPCA/functions')
addpath('JisstPCA/simulations/functions')


p = 80; q_vec = [50 80]; N_vec = [20 40];
rx = [3 2]; ry = rx;
K = 2;

u = cell(K, 1); V = cell(K, 1); W = cell(K, 1);
filename = "JisstPCA/simulations/ComparativeStudies/NetworkSim/results.csv";
cluster_filename = "JisstPCA/simulations/ComparativeStudies/NetworkSim/results_cluster.csv";
truecluster_path = "JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_";
for seed = 1 : 20
    for q = q_vec
        % probability tensor construction
        px1 = table2array(readtable(sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_px1_q%d.csv", q)));
        px2 = table2array(readtable(sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_px2_q%d.csv", q)));
        py1 = table2array(readtable(sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_py1_q%d.csv", q)));
        py2 = table2array(readtable(sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/FactorsRandom/setting4_py2_q%d.csv", q)));
        thr = 0.75;
        for N = N_vec
            rng(seed)
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
            [V{1},Dx{1},~] =svds(px_1, rx(1) - 1); [V{2},Dx{2},~] =svds(px_2, rx(2) - 1);
            [W{1},Dy{1},~] =svds(py_1, ry(1) - 1); [W{2},Dy{2},~] =svds(py_2, ry(2) - 1);
            Network_x{1} = px_1 * norm(u1); Network_x{2} = px_2 * norm(u2);
            Network_y{1} = py_1 * norm(u1); Network_y{2} = py_2 * norm(u2);
            u{1} = u1/norm(u1); u{2} = u2/norm(u2);
            
            vnames = {'q', 'N', 'seed', 'method', 'factor name', 'factor number', 'err'};
            results = array2table(zeros(0,7), 'VariableNames',vnames);
            vnames = {'q', 'N', 'seed', 'method', 'factor name', 'RI', 'ARI'};
            cluster_results = array2table(zeros(0,7), 'VariableNames',vnames);

            % data tensors
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
            VW_results = array2table(zeros(0, 6), 'VariableNames', vnames); VW_weighted_results =  VW_results;
            vnames = {'method', 'FactorNumber', 'RowInd', 'Val'};
            u_results = array2table(zeros(0, 4), 'VariableNames', vnames);

            VW_results = AddStructuredValVW(VW_results, V, W, "truth"); 
            u_results = AddStructuredValU(u_results, u, "truth");

            
            %subtraction deflation & partial projection (sample) deflation
            for deflation = [0 2]
                %Jisst-PCA, oracle ranks
                tic
                [u_est, V_est, W_est, d_est] = JisstPCA(X_obs, Y_obs, 2, 'rx', rx - 1, 'ry', ry - 1,...
                    'deflation', deflation, 'max_iter', 10);
                toc 
                Dx_est = cell(K, 1); Dy_est = cell(K, 1);
                for k = 1 : K
                    Dx_est{k} = ones(size(V_est{k}, 2)) * d_est(1, k); Dy_est{k} = ones(size(W_est{k}, 2)) * d_est(2, k);
                end
                results = add_results_structured(results, q, N, seed, sprintf("JisstPCA%d (oracle)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                VW_results = AddStructuredValVW(VW_results, V_est, W_est, sprintf("JisstPCA%d (oracle)", deflation));
                u_results = AddStructuredValU(u_results, u_est, sprintf("JisstPCA%d (oracle)", deflation)); 
                [RI, ARI] = find_cluster_ARI(u_est, V_est, W_est, Dx_est, Dy_est, q, u, [rx;ry], truecluster_path);
                cluster_results = add_results_cluster(cluster_results, q, N, seed, sprintf("JisstPCA%d (oracle)", deflation), RI, ARI);
            
                %Jisst-PCA, BIC ranks
                tic
                [u_est, V_est, W_est, d_est] = JisstPCA(X_obs, Y_obs, 2, 'deflation', deflation, 'max_iter', 10);
                toc 
                Dx_est = cell(K, 1); Dy_est = cell(K, 1);
                for k = 1 : K
                    Dx_est{k} = ones(size(V_est{k}, 2)) * d_est(1, k); Dy_est{k} = ones(size(W_est{k}, 2)) * d_est(2, k);
                end
                results = add_results_structured(results, q, N, seed, sprintf("JisstPCA%d (BIC)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                VW_results = AddStructuredValVW(VW_results, V_est, W_est, sprintf("JisstPCA%d (BIC)", deflation));
                u_results = AddStructuredValU(u_results, u_est, sprintf("JisstPCA%d (BIC)", deflation)); 
                [RI, ARI] = find_cluster_ARI(u_est, V_est, W_est, Dx_est, Dy_est, q, u, [rx;ry], truecluster_path);
                cluster_results = add_results_cluster(cluster_results, q, N, seed, sprintf("JisstPCA%d (BIC)", deflation), RI, ARI);

                %G-JisstPCA, oracle ranks
                tic
                [u_est, V_est, W_est, Dx_est, Dy_est] = dJisstPCA(X_obs, Y_obs, 2, 'rx', rx - 1, 'ry', ry - 1, ...
                    'deflation', deflation, 'max_iter', 10);
                toc
                results = add_results_structured(results, q, N, seed, sprintf("G-JisstPCA%d (oracle)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                VW_results = AddStructuredValVW(VW_results, V_est, W_est, sprintf("G-JisstPCA%d (oracle)", deflation));
                u_results = AddStructuredValU(u_results, u_est, sprintf("G-JisstPCA%d (oracle)", deflation)); 
                [RI, ARI] = find_cluster_ARI(u_est, V_est, W_est, Dx_est, Dy_est, q, u, [rx;ry], truecluster_path);
                cluster_results = add_results_cluster(cluster_results, q, N, seed, sprintf("G-JisstPCA%d (oracle)", deflation), RI, ARI);

                %G-JisstPCA, BIC ranks
                tic
                [u_est, V_est, W_est, Dx_est, Dy_est] = dJisstPCA(X_obs, Y_obs, 2, 'deflation', deflation, 'max_iter', 10);
                toc
                results = add_results_structured(results, q, N, seed, sprintf("G-JisstPCA%d (BIC)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                VW_results = AddStructuredValVW(VW_results, V_est, W_est, sprintf("G-JisstPCA%d (BIC)", deflation));
                u_results = AddStructuredValU(u_results, u_est, sprintf("G-JisstPCA%d (BIC)", deflation)); 
                [RI, ARI] = find_cluster_ARI(u_est, V_est, W_est, Dx_est, Dy_est, q, u, [rx;ry], truecluster_path);
                cluster_results = add_results_cluster(cluster_results, q, N, seed, sprintf("G-JisstPCA%d (BIC)", deflation), RI, ARI);
            end

            % Method 3: iHOSVD
            [u_est_mat, V_est_mat, W_est_mat, ~, ~] = iHOSVD(X_obs, Y_obs, sum(rx - 1), sum(ry - 1), 2);
            u_est = cell(K, 1); V_est = cell(K, 1); W_est = cell(K, 1); start_V = 1; start_W = 1;
            for k = 1 : K
                u_est{k} = u_est_mat(:, k);
                V_est{k} = V_est_mat(:, start_V : (start_V + rx(k) - 2));
                W_est{k} = W_est_mat(:, start_W : (start_W + ry(k) - 2));
                start_V = start_V + rx(k) - 1; start_W = start_W + ry(k) - 1;
            end
            results = add_results_structured(results, q, N, seed, "iHOSVD", find_dist(u_est, V_est, W_est, u, V, W));
            VW_results = AddStructuredValVW(VW_results, V_est, W_est, "iHOSVD");
            u_results = AddStructuredValU(u_results, u_est, "iHOSVD"); 
            for k = 1 : K
                Dx_est{k} = ones(size(V_est{k}, 2)); Dy_est{k} = ones(size(W_est{k}, 2));
            end
            [RI, ARI] = find_cluster_ARI(u_est, V_est, W_est, Dx_est, Dy_est, q, u, [rx;ry], truecluster_path);
            cluster_results = add_results_cluster(cluster_results, q, N, seed, "iHOSVD", RI, ARI);


            % Method 4: iHOOI
            [u_est_mat, V_est_mat, W_est_mat, ~, ~] = iHOOI(X_obs, Y_obs, sum(rx - 1), sum(ry - 1), 10, 2);
            u_est = cell(K, 1); V_est = cell(K, 1); W_est = cell(K, 1); start_V = 1; start_W = 1;
            for k = 1 : K
                u_est{k} = u_est_mat(:, k);
                V_est{k} = V_est_mat(:, start_V : (start_V + rx(k) - 2));
                W_est{k} = W_est_mat(:, start_W : (start_W + ry(k) - 2));
                start_V = start_V + rx(k) - 1; start_W = start_W + ry(k) - 1;
            end
            results = add_results_structured(results, q, N, seed, "iHOOI", find_dist(u_est, V_est, W_est, u, V, W));
            VW_results = AddStructuredValVW(VW_results, V_est, W_est, "iHOOI");
            u_results = AddStructuredValU(u_results, u_est, "iHOOI"); 
            for k = 1 : K
                Dx_est{k} = ones(size(V_est{k}, 2)); Dy_est{k} = ones(size(W_est{k}, 2));
            end
            [RI, ARI] = find_cluster_ARI(u_est, V_est, W_est, Dx_est, Dy_est, q, u, [rx;ry], truecluster_path);
            cluster_results = add_results_cluster(cluster_results, q, N, seed, "iHOOI", RI, ARI);

            VW_weighted_filename = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/results/factorVW_weighted_q%d_N%d_seed%d.csv", q, N, seed);
            VW_filename = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/results/factorVW_q%d_N%d_seed%d.csv", q, N, seed);
            u_filename = sprintf("JisstPCA/simulations/ComparativeStudies/NetworkSim/results/factorU_q%d_N%d_seed%d.csv", q, N, seed);
            
            if isfile(filename)
               writetable(results, filename, 'WriteMode', 'Append', 'WriteVariableNames', false);
               writetable(VW_results, VW_filename, 'WriteMode', 'Append', 'WriteVariableNames', false);
               writetable(u_results, u_filename, 'WriteMode', 'Append', 'WriteVariableNames', false);
               writetable(cluster_results, cluster_filename, 'WriteMode', 'Append', 'WriteVariableNames', false);
            else
               writetable(results, filename);
               writetable(VW_results, VW_filename);
               writetable(u_results, u_filename);
               writetable(cluster_results, cluster_filename);
            end
            sprintf("finished q = %d, N = %d, seed = %d", q, N, seed)
        end
    end
end

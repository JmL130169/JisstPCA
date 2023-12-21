%% comparative study, unstructured data, scalar d setting
% Software dependency: tensor toolbox (https://www.tensortoolbox.org)
addpath('JisstPCA/functions')

scenario_list = 1 : 4; %4 scenarios of dimensions
orthogonal_list = [true false];
SNR_list = [4 5 6.25 8.5 12.5 25];
p_vec = [50 150 150 150]; q_vec = [50 150 50 50]; N_vec = [200 50 200 50];
nrep = 10;

rx = [3 2]; ry = rx;
K = 2;
scenario = 4;

for seed = 1 : 10
    for setting = 1 : 2
        filename = sprintf("JisstPCA/simulations/ComparativeStudies/UnstructuredSim/scalar_d_res%d.csv", setting);
        for SNR_k = 1 : length(SNR_list)
            for orth_k = 1 : 2
                vnames = {'scenario', 'orthogonal', 'SNR', 'seed', 'method', 'factor name', 'factor number', 'err'};
                results = array2table(zeros(0,8), 'VariableNames',vnames);

                SNR = SNR_list(SNR_k);
                p = p_vec(scenario);
                q = q_vec(scenario);
                N = N_vec(scenario);
                if setting == 1
                    dx0 = SNR * (sqrt(p) + sqrt(N)); dx = [dx0 dx0/2]; 
                    dy0 = SNR * (sqrt(q) + sqrt(N)); dy = [dy0 dy0/2];
                else
                    dx0 = SNR * (sqrt(p) + sqrt(N)); dx = [dx0 dx0]; 
                    dy0 = SNR * (sqrt(q) + sqrt(N)); dy = [dy0 dy0*9/10];
                end
                Dx = cell(K, 1); Dy = cell(K, 1);
                for k = 1 : K
                    Dx{k} = dx(k) * diag(ones(rx(k),1));
                    Dy{k} = dy(k) * diag(ones(ry(k),1));
                end

                %model generation
                %rng(seed)
                rng(1234)
                [X, Y, u, V, W] = gen_sstTensors(p, q, N, rx, ry, K, orthogonal_list(orth_k), Dx, Dy);
                %data generation
                rng(seed)
                epi_X = wigner_sst(p, N);
                epi_Y = wigner_sst(q, N);
                X_obs = tensor(double(X) + epi_X);
                Y_obs = tensor(double(Y) + epi_Y);

            
                for deflation = 0 : 1
                    %Jisst-PCA, BIC & oracle ranks, orthogonal & subtraction deflation
                    tic
                    [u_est, V_est, W_est, ~] = JisstPCA(X_obs, Y_obs, 2, 'deflation', deflation, 'max_iter', 10);
                    toc 
                    results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                    sprintf("JisstPCA%d (BIC)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                
                    V_est_top = cell(K, 1); W_est_top = cell(K, 1);
                    for k = 1 : K
                        rx_tmp = min(size(V_est{k}, 2), rx(k));
                        ry_tmp = min(size(W_est{k}, 2), ry(k));
                        V_est_top{k} = V_est{k}(:, 1 : rx_tmp); W_est_top{k} = W_est{k}(:, 1 : ry_tmp);
                    end
                    results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                        sprintf("JisstPCA%d (BIC top)", deflation), find_dist(u_est, V_est_top, W_est_top, u, V, W));
               
                    if max(find_ranks(V_est, 2) ~= rx) || max(find_ranks(W_est, 2) ~= ry)
                        [u_est, V_est, W_est, ~] = JisstPCA(X_obs, Y_obs, 2, 'rx', rx, 'ry', ry, 'deflation', deflation, 'max_iter', 10);          
                    end
                    results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                        sprintf("JisstPCA%d (oracle)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                    sprintf("Finished scenario %d, orthogonal = %d, SNR index = %d, seed = %d, deflation = %d, JisstPCA", scenario, orthogonal_list(orth_k), SNR_k, seed, deflation)
                
        
                    %G-JisstPCA
                    tic
                    [u_est, V_est, W_est, Dx_est, Dy_est] = dJisstPCA(X_obs, Y_obs, 2, 'deflation', deflation, 'max_iter', 10);
                    toc
                    results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                       sprintf("G-JisstPCA%d (BIC)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                    V_est_top = cell(K, 1); W_est_top = cell(K, 1);
                    for k = 1 : K
                        rx_tmp = min(size(V_est{k}, 2), rx(k));
                        ry_tmp = min(size(W_est{k}, 2), ry(k));
                        V_est_top{k} = V_est{k}(:, 1 : rx_tmp); W_est_top{k} = W_est{k}(:, 1 : ry_tmp);
                    end
                    results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                        sprintf("G-JisstPCA%d (BIC top)", deflation), find_dist(u_est, V_est_top, W_est_top, u, V, W));
               
                    if max(find_ranks(V_est, 2) ~= rx) || max(find_ranks(W_est, 2) ~= ry)
                        [u_est, V_est, W_est, ~] = dJisstPCA(X_obs, Y_obs, 2, 'rx', rx, 'ry', ry, 'deflation', deflation, 'max_iter', 10);
                    end
                    results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                        sprintf("G-JisstPCA%d (oracle)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
                    sprintf("Finished scenario %d, orthogonal = %d, SNR index = %d, seed = %d, deflation = %d, G-JisstPCA", scenario, orthogonal_list(orth_k), SNR_k, seed, deflation)

                end

                % Method 3: iHOSVD
                [u_est, V_est, W_est, ~, ~] = iHOSVD(X_obs, Y_obs, sum(rx), sum(ry), 2);
                results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                    "iHOSVD", find_dist_mat(u_est, V_est, W_est, u, V, W));
        

                % Method 4: iHOOI
                [u_est, V_est, W_est, ~, ~] = iHOOI(X_obs, Y_obs, sum(rx), sum(ry), 10, 2);
                results = add_results_unstructured(results, scenario, orthogonal_list(orth_k), SNR, seed,...
                    "iHOOI", find_dist_mat(u_est, V_est, W_est, u, V, W));
    
                if isfile(filename)
                    writetable(results, filename,'WriteMode','Append','WriteVariableNames',false);
                else
                    writetable(results, filename);
                end
            end
        end
    end
end

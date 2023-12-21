%% comparative study, structured data with diagonal D (eigenvalues more different than previous code)
% Software dependency: tensor toolbox (https://www.tensortoolbox.org)
addpath('JisstPCA/functions')
filename = "JisstPCA/simulations/ComparativeStudies/GeneralizedModel/diagnoal_d_res4.csv";


SNR_list = [4 5 6.25 8.5 12.5 25];
nrep = 10;
p = 150; q = 50; N = 50; 
rx = [3 2]; ry = rx;
K = 2;
structure_list = ["block" "star"]; sample_clusters = 3;

for seed = 1 : 10
    for SNR_k = 1 : 6
        vnames = {'diagonal', 'SNR', 'seed', 'method', 'factor name', 'factor number', 'err'};
        results = array2table(zeros(0,7), 'VariableNames',vnames);

        SNR = SNR_list(SNR_k);
        dx0 = SNR * (sqrt(p) + sqrt(N)); dx = [dx0 0.5 * dx0]; 
        dy0 = SNR * (sqrt(q) + sqrt(N)); dy = [dy0 0.5 * dy0];
        Dx = cell(K, 1); Dy = cell(K, 1);
        Dx{1} = dx(1) * diag([3.2 2 1.2]); Dx{2} = dx(2) * diag([2 1.6]);
        Dy{1} = dy(1) * diag([3.2 2 1.2]); Dy{2} = dy(2) * diag([2 1.6]);

        rng(1234)
        [X, Y, u, V, W] = gen_StructTensors(p, q, N, rx, ry, K, structure_list, Dx, Dy, sample_clusters);
        %plot(u{1}, u{2}, 'o')
        %colormap('parula')
        %imagesc(abs(V{1}*V{1}'))
        %colorbar
        
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
            results = add_results_structured(results, diagonal, SNR, seed,...
            sprintf("JisstPCA%d (BIC)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
            
            V_est_top = cell(K, 1); W_est_top = cell(K, 1);
            for k = 1 : K
                rx_tmp = min(size(V_est{k}, 2), rx(k));
                ry_tmp = min(size(W_est{k}, 2), ry(k));
                V_est_top{k} = V_est{k}(:, 1 : rx_tmp); W_est_top{k} = W_est{k}(:, 1 : ry_tmp);
            end
            results = add_results_structured(results, diagonal, SNR, seed,...
                    sprintf("JisstPCA%d (BIC top)", deflation), find_dist(u_est, V_est_top, W_est_top, u, V, W));
               
            if max(find_ranks(V_est, 2) ~= rx) || max(find_ranks(W_est, 2) ~= ry)
                [u_est, V_est, W_est, ~] = JisstPCA(X_obs, Y_obs, 2, 'rx', rx, 'ry', ry, 'deflation', deflation, 'max_iter', 10);
            end
            results = add_results_structured(results, diagonal, SNR, seed,...
                sprintf("JisstPCA%d (oracle)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
            sprintf("Finished diagonal = %d, SNR index = %d, seed = %d, deflation = %d, JisstPCA", diagonal, SNR_k, seed, deflation)
            
            
            %G-JisstPCA
            tic
            [u_est, V_est, W_est, Dx_est, Dy_est] = dJisstPCA(X_obs, Y_obs, 2, 'deflation', deflation, 'max_iter', 10);
            toc
            results = add_results_structured(results, diagonal, SNR, seed,...
                sprintf("G-JisstPCA%d (BIC)", deflation), find_dist(u_est, V_est, W_est, u, V, W));

            V_est_top = cell(K, 1); W_est_top = cell(K, 1);
            for k = 1 : K
                rx_tmp = min(size(V_est{k}, 2), rx(k));
                ry_tmp = min(size(W_est{k}, 2), ry(k));
                V_est_top{k} = V_est{k}(:, 1 : rx_tmp); W_est_top{k} = W_est{k}(:, 1 : ry_tmp);
            end
            results = add_results_structured(results, diagonal, SNR, seed,...
                    sprintf("G-JisstPCA%d (BIC top)", deflation), find_dist(u_est, V_est_top, W_est_top, u, V, W));
                             
            if max(find_ranks(V_est, 2) ~= rx) || max(find_ranks(W_est, 2) ~= ry)
                [u_est, V_est, W_est, ~] = dJisstPCA(X_obs, Y_obs, 2, 'rx', rx, 'ry', ry, 'deflation', deflation, 'max_iter', 10);
            end
            results = add_results_structured(results, diagonal, SNR, seed,...
                sprintf("G-JisstPCA%d (oracle)", deflation), find_dist(u_est, V_est, W_est, u, V, W));
            sprintf("Finished diagonal = %d, SNR index = %d, seed = %d, deflation = %d, G-JisstPCA", diagonal, SNR_k, seed, deflation)
        end

        % Method 3: iHOSVD
        [u_est, V_est, W_est, ~, ~] = iHOSVD(X_obs, Y_obs, sum(rx), sum(ry), 2);
        results = add_results_structured(results, diagonal, SNR, seed,...
            "iHOSVD", find_dist_mat(u_est, V_est, W_est, u, V, W));
    

        % Method 4: iHOOI
        [u_est, V_est, W_est, ~, ~] = iHOOI(X_obs, Y_obs, sum(rx), sum(ry), 10, 2);
        results = add_results_structured(results, diagonal, SNR, seed,...
            "iHOOI", find_dist_mat(u_est, V_est, W_est, u, V, W));
    
        if isfile(filename)
            writetable(results, filename,'WriteMode','Append','WriteVariableNames',false);
        else
            writetable(results, filename);
        end
    end
end

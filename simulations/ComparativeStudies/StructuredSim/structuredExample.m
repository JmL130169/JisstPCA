%% comparative study, structured data example
% Software dependency: tensor toolbox (https://www.tensortoolbox.org)
addpath('JisstPCA/functions')


SNR = 3;
p = 50; q = 50; N = 200; 
rx = [3 2]; ry = rx;
K = 2;
structure_list = ["block" "star"]; sample_clusters = 3;


seed = 1;
vnames = {'diagonal', 'SNR', 'seed', 'method', 'factor name', 'factor number', 'err'};
results = array2table(zeros(0,7), 'VariableNames',vnames);

dx0 = SNR * (sqrt(p) + sqrt(N)); dx = [dx0 0.5 * dx0]; 
dy0 = SNR * (sqrt(q) + sqrt(N)); dy = [dy0 0.5 * dy0];
Dx = cell(K, 1); Dy = cell(K, 1);

for k = 1 : K
    Dx{k} = dx(k) * diag(ones(rx(k),1));
    Dy{k} = dy(k) * diag(ones(ry(k),1));
end

u = cell(K, 1); V = cell(K, 1); W = cell(K, 1);
u{1} = table2array(readtable("JisstPCA/simulations/ComparativeStudies/StructuredSims/u1.csv"));
u{2} = table2array(readtable("JisstPCA/simulations/ComparativeStudies/StructuredSims/u2.csv"));
V{1} = table2array(readtable("JisstPCA/simulations/ComparativeStudies/StructuredSims/V1.csv"));
V{2} = table2array(readtable("JisstPCA/simulations/ComparativeStudies/StructuredSims/V2.csv"));
W{1} = table2array(readtable("JisstPCA/simulations/ComparativeStudies/StructuredSims/W1.csv"));
W{2} = table2array(readtable("JisstPCA/simulations/ComparativeStudies/StructuredSims/W2.csv"));
p = size(V{1}, 1); q = size(W{1}, 1); N = length(u{1});
    X = tensor(zeros(p, p, N)); Y = tensor(zeros(q, q, N)); 
    for k = 1 : K
        X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(u{k})));
        Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(u{k})));
    end
vnames = {'method', 'FactorName', 'FactorNumber', 'RowInd', 'ColInd', 'Val'};
VW_results = array2table(zeros(0, 6), 'VariableNames', vnames);
vnames = {'method', 'FactorNumber', 'RowInd', 'Val'};
u_results = array2table(zeros(0, 4), 'VariableNames', vnames);

VW_results = AddStructuredValVW(VW_results, V, W, "truth");
u_results = AddStructuredValU(u_results, u, "truth");
rng(seed)
epi_X = wigner_sst(p, N);
epi_Y = wigner_sst(q, N);
X_obs = tensor(double(X) + epi_X);
Y_obs = tensor(double(Y) + epi_Y);

deflation = 0;
%Jisst-PCA, BIC & oracle ranks, orthogonal & subtraction deflation
tic
[u_est, V_est, W_est, ~] = JisstPCA(X_obs, Y_obs, 2, 'deflation', deflation, 'max_iter', 10);
toc 
VW_results = AddStructuredValVW(VW_results, V_est, W_est, "JisstPCA");
u_results = AddStructuredValU(u_results, u_est, "JisstPCA");        
        
%G-JisstPCA
tic
[u_est, V_est, W_est, Dx_est, Dy_est] = dJisstPCA(X_obs, Y_obs, 2, 'deflation', deflation, 'max_iter', 10);
toc
VW_results = AddStructuredValVW(VW_results, V_est, W_est, 'G-JisstPCA');
u_results = AddStructuredValU(u_results, u_est, 'G-JisstPCA'); 


% Method 3: iHOSVD
[u_est_mat, V_est_mat, W_est_mat, ~, ~] = iHOSVD(X_obs, Y_obs, sum(rx), sum(ry), 2);
u_est = cell(K, 1); V_est = cell(K, 1); W_est = cell(K, 1); start_V = 1; start_W = 1;
for k = 1 : K
    u_est{k} = u_est_mat(:, k);
    V_est{k} = V_est_mat(:, start_V : (start_V + rx(k) - 1));
    W_est{k} = W_est_mat(:, start_W : (start_W + ry(k) - 1));
    start_V = start_V + rx(k); start_W = start_W + ry(k);
end
VW_results = AddStructuredValVW(VW_results, V_est, W_est, "iHOSVD");
u_results = AddStructuredValU(u_results, u_est, "iHOSVD"); 

%[u_est, V_est, W_est, ~] = dJisstPCA(X_obs, Y_obs, 2, 'rx', rx, 'ry', ry, 'deflation', deflation, 'max_iter', 10);
%find_dist(u_est, V_est, W_est, u, V, W)
%colormap('parula')
%imagesc(abs(V_est{2}*V_est{2}'))

% Method 4: iHOOI
[u_est, V_est, W_est, ~, ~] = iHOOI(X_obs, Y_obs, sum(rx), sum(ry), 10, 2);
u_est = cell(K, 1); V_est = cell(K, 1); W_est = cell(K, 1); start_V = 1; start_W = 1;
for k = 1 : K
    u_est{k} = u_est_mat(:, k);
    V_est{k} = V_est_mat(:, start_V : (start_V + rx(k) - 1));
    W_est{k} = W_est_mat(:, start_W : (start_W + ry(k) - 1));
    start_V = start_V + rx(k); start_W = start_W + ry(k);
end
%find_dist_mat(u_est, V_est, W_est, u, V, W)
%colormap('parula')
%imagesc(abs(V_est{2} * V_est{2}'))
VW_results = AddStructuredValVW(VW_results, V_est, W_est, "iHOOI");
u_results = AddStructuredValU(u_results, u_est, "iHOOI"); 
        

VW_filename = "JisstPCA/simulations/ComparativeStudies/StructuredSims/factorVW.csv";
u_filename = "JisstPCA/simulations/ComparativeStudies/StructuredSims/factorU.csv";
writetable(VW_results, VW_filename);
writetable(u_results, u_filename);

%% comparative study, network data
mkdir('/Users/jiamingliu/Desktop/tensor_toolbox')   
addpath('/Users/jiamingliu/Desktop/tensor_toolbox')
mkdir('/Users/jiamingliu/Desktop/isst')   
addpath('/Users/jiamingliu/Desktop/isst')

%% case1: p = q = 150, N = 200
p = 150;
q = 150;
N = 200;
r = [3 2];
tol = 0.0001;
max_iter = 100;
d = [1 1]; % contribution of each SBM

% SBM setup
px1 = readtable("/Users/jiamingliu/Desktop/section 4/network_data/case1/px1.csv");
px1 = double(table2array(px1(:, 2:(p+1))));
px2 = readtable("/Users/jiamingliu/Desktop/section 4/network_data/case1/px2.csv");
px2 = double(table2array(px2(:, 2:(p+1))));
py1 = readtable("/Users/jiamingliu/Desktop/section 4/network_data/case1/py1.csv");
py1 = double(table2array(py1(:, 2:(p+1))));
py2 = readtable("/Users/jiamingliu/Desktop/section 4/network_data/case1/py2.csv");
py2 = double(table2array(py2(:, 2:(p+1))));

% prepocessing (true probability matrix)
n_x1 = ones(p, 1);
n_y1 = ones(q, 1);
px_1 = (eye(p)-n_x1*n_x1'/p)*px1*(eye(p)-n_x1*n_x1'/p);
px_2 = (eye(p)-n_x1*n_x1'/p)*px2*(eye(p)-n_x1*n_x1'/p);
py_1 = (eye(q)-n_y1*n_y1'/q)*py1*(eye(q)-n_y1*n_y1'/q);
py_2 = (eye(q)-n_y1*n_y1'/q)*py2*(eye(q)-n_y1*n_y1'/q);

r1 = r - 1;

[a1, b1, ~] = svd(px_1);
b1 = diag(b1);
Dx1 = diag(b1(1:r1(1))); % true Dx1
V1 = a1(:, 1:r1(1)); % true V1
[a2, b2, ~] = svd(px_2);
b2 = diag(b2);
Dx2 = diag(b2(1:r1(2))); % true Dy2
V2 = a2(:, 1:r1(2)); % true V2

[a3, b3, ~] = svd(py_1);
b3 = diag(b3);
Dy1 = diag(b3(1:r1(1))); % true Dy1
W1 = a3(:, 1:r1(1)); % true W1
[a4, b4, ~] = svd(py_2);
b4 = diag(b4);
Dy2 = diag(b4(1:r1(2))); % true Dy2
W2 = a4(:, 1:r1(2)); % true W2

% tensor construction
rng(2023)
thr = 0.75;
u1 = zeros(N, 1);
u1_norm = u1/norm(u1);
for i = 1:N
    ur = rand();
    if (ur < thr)
        u1(i) = 1;
    else
        u1(i) = 0;
    end
end
u2 = 1 - u1;
u2_norm = u2/norm(u2);

% true probability tensors
X1 = squeeze(ttt(tensor(px1), tensor(u1)));
X2 = squeeze(ttt(tensor(px2), tensor(u2)));
Y1 = squeeze(ttt(tensor(py1), tensor(u1)));
Y2 = squeeze(ttt(tensor(py2), tensor(u2)));
X_t = d(1)*X1 + d(2)*X2;
Y_t = d(1)*Y1 + d(2)*Y2;

rep = 10;
err_jisst = cell(3, 1);
err_ihosvd = cell(3, 1);
err_ihooi = cell(3, 1);
for i = 1:3
    err_jisst{i} = zeros(rep, 2);
    err_ihosvd{i} = zeros(rep, 2);
    err_ihooi{i} = zeros(rep, 2);
end

% generate rep = 10 adjacency tensors and estimate its error
Ax = cell(rep, 1);
Ay = cell(rep, 1);
Ax_shf = cell(rep, 1);
Ay_shf = cell(rep, 1);
u0 = cell(rep, 1);

for l = 1:rep

    rng(l)
    Ax{l} = zeros(size(X_t));
    Ay{l} = zeros(size(Y_t));
    for k = 1:N
        for j = 1:p
            for i = 1:(j-1)
                u_rand = rand();
                if(X_t(i, j, k) <= u_rand)
                    Ax{l}(i, j, k) = 0;
                else
                    Ax{l}(i, j, k) = 1;
                end
                Ax{l}(j, i, k) = Ax{l}(i, j, k);
            end
        end
    end
    for k = 1:N
        for j = 1:q
            for i = 1:(j-1)
                u_rand = rand();
                if(Y_t(i, j, k) <= u_rand)
                    Ay{l}(i, j, k) = 0;
                else
                    Ay{l}(i, j, k) = 1;
                end
                Ay{l}(j, i, k) = Ay{l}(i, j, k);
            end
        end
    end

    Ax_shf{l} = zeros(size(Ax{l}));
    Ay_shf{l} = zeros(size(Ay{l}));
    for k = 1:N
        Ax_shf{l}(:, :, k) = (eye(p)-n_x1*n_x1'/p)*Ax{l}(:, :, k)*(eye(p)-n_x1*n_x1'/p);
        Ay_shf{l}(:, :, k) = (eye(q)-n_y1*n_y1'/q)*Ay{l}(:, :, k)*(eye(q)-n_y1*n_y1'/q);
    end
    Ax_shf{l} = tensor(Ax_shf{l});
    Ay_shf{l} = tensor(Ay_shf{l});

    % initialization
    M = [double(tenmat(Ax_shf{l}, 3, [1, 2])), double(tenmat(Ay_shf{l}, 3, [1, 2]))];
    [a, b, ~] = svd(M);
    [~, ind] = sort(diag(b));
    a = a(:, ind);
    a = fliplr(a);
    hat_u = a(:, 1);
    u0{l} = hat_u; % spectral initialization for adjacnecy tensor

    lam = norm(Ax_shf{l})/(norm(Ax_shf{l}) + norm(Ay_shf{l}));
    lam_adj = [lam lam];

    % estimation using different methods and report the errors
    % JisstPCA_diag
    [u1_est, V1_est, W1_est, ~, ~] = mjisst_PCA_sub(Ax_shf{l}, Ay_shf{l}, u0{l}, r1, r1, lam_adj, tol, max_iter);
    
    err_jisst{1}(l, 1) = abs(sin_do(u1_est{2}, u1_norm));
    err_jisst{1}(l, 2) = abs(sin_do(u1_est{3}, u2_norm));
    err_jisst{2}(l, 1) = abs(sin_do(V1_est{2}, V1));
    err_jisst{2}(l, 2) = abs(sin_do(V1_est{3}, V2));
    err_jisst{3}(l, 1) = abs(sin_do(W1_est{2}, W1));
    err_jisst{3}(l, 2) = abs(sin_do(W1_est{3}, W2));

    % iHOSVD
    [u2_est, V2_est, W2_est, ~, ~] = iHOSVD(Ax_shf{l}, Ay_shf{l}, sum(r1), sum(r1), 2);

    err_ihosvd{1}(l, 1) = abs(sin_do(u1_norm, u2_est(:, 1)));
    err_ihosvd{1}(l, 2) = abs(sin_do(u2_norm, u2_est(:, 2)));
    err_ihosvd{2}(l, 1) = abs(sin_do(V1, V2_est(:, 1:2)));
    err_ihosvd{2}(l, 2) = abs(sin_do(V2, V2_est(:, 3)));
    err_ihosvd{3}(l, 1) = abs(sin_do(W1, W2_est(:, 1:2)));
    err_ihosvd{3}(l, 2) = abs(sin_do(W2, W2_est(:, 3)));

    % iHOOI
    [u3_est, V3_est, W3_est, ~, ~] = iHOOI(Ax_shf{l}, Ay_shf{l}, sum(r1), sum(r1), max_iter, 2);

    err_ihooi{1}(l, 1) = abs(sin_do(u1_norm, u3_est(:, 1)));
    err_ihooi{1}(l, 2) = abs(sin_do(u2_norm, u3_est(:, 2)));
    err_ihooi{2}(l, 1) = abs(sin_do(V1, V3_est(:, 1:2)));
    err_ihooi{2}(l, 2) = abs(sin_do(V2, V3_est(:, 3)));
    err_ihooi{3}(l, 1) = abs(sin_do(W1, W3_est(:, 1:2)));
    err_ihooi{3}(l, 2) = abs(sin_do(W2, W3_est(:, 3)));
end

% export the data
writematrix(err_jisst{1}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_jisst_u.csv")
writematrix(err_jisst{2}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_jisst_V.csv")
writematrix(err_jisst{3}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_jisst_W.csv")

writematrix(err_ihosvd{1}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_ihosvd_u.csv")
writematrix(err_ihosvd{2}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_ihosvd_V.csv")
writematrix(err_ihosvd{3}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_ihosvd_W.csv")

writematrix(err_ihooi{1}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_ihooi_u.csv")
writematrix(err_ihooi{2}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_ihooi_V.csv")
writematrix(err_ihooi{3}, "/Users/jiamingliu/Desktop/section 4/network_data/case1/err_ihooi_W.csv")











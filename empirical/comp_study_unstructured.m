%% comparative study, unstructured data
mkdir('/Users/jl259/Desktop/tensor_toolbox')   
addpath('/Users/jl259/Desktop/tensor_toolbox')
mkdir('/Users/jl259/Desktop/isst')   
addpath('/Users/jl259/Desktop/isst')

%% scenario 1: p = q = 50, N = 200
p = 50;
q = 50;
N = 200;
r = [3 2];
dx = [2*p p];
dy = [2*q q];
tol = 0.0001;
max_iter = 100;

rep = 10;
eta = 0.4:0.4:2; % SNR level
sz = size(eta, 2);

err_oracle = cell(3, 2);
err_bic = cell(3, 2);
err_ihosvd = cell(3, 2);
err_ihooi = cell(3, 2);
for i = 1:3
    for j = 1:2
        err_oracle{i, j} = zeros(rep, sz);
        err_bic{i, j} = zeros(rep, sz);
        err_ihosvd{i, j} = zeros(rep, sz);
        err_ihooi{i, j} = zeros(rep, sz);
    end
end

rng("default")
V1 = orth(normrnd(0, 1, [p, r(1)]));
V2 = orth(normrnd(0, 1, [p, r(2)]));
W1 = orth(normrnd(0, 1, [q, r(1)]));
W2 = orth(normrnd(0, 1, [q, r(2)]));
u = normrnd(0, 1, [N, 3]);
u1 = u(:, 1)/norm(u(:, 1));
u2 = u(:, 2)/norm(u(:, 2));
X1 = squeeze(ttt(tensor(V1*V1'), tensor(u1)));
X2 = squeeze(ttt(tensor(V2*V2'), tensor(u2)));
Y1 = squeeze(ttt(tensor(W1*W1'), tensor(u1)));
Y2 = squeeze(ttt(tensor(W2*W2'), tensor(u2)));
X_t = dx(1)*X1 + dx(2)*X2;
Y_t = dy(1)*Y1 + dy(2)*Y2;

for i = 1:rep
    
    % set noise and perturbed tensor
    epi_X = zeros(size(X_t));
    epi_Y = zeros(size(Y_t));

    rng(i)
    for k = 1:N
        epi_X(:, :, k) = wigner(p, 1);
        epi_Y(:, :, k) = wigner(q, 1);
    end
    
    X = cell(sz, 1); 
    Y = cell(sz, 1);
    for j = 1:sz
        X{j} = tensor(eta(j)*(norm(tensor(epi_X))/norm(X_t))*X_t + epi_X);
        Y{j} = tensor(eta(j)*(norm(tensor(epi_Y))/norm(Y_t))*Y_t + epi_Y);
    end

    for j = 1:sz

        % initialization
        M = [double(tenmat(X{j}, 3, [1, 2])), double(tenmat(Y{j}, 3, [1, 2]))];
        [a, b, ~] = svd(M);
        [~, ind] = sort(diag(b));
        a = a(:, ind);
        a = fliplr(a);
        hat_u = a(:, 1);
        u0 = hat_u; % spectral initialization

        % scaler lambda
        lam = norm(X{j})/(norm(X{j}) + norm(Y{j}));
        lambda = [lam lam];

        % sum/cumsum of r
        r_sum = sum(r);
        cr = cumsum(r);

        % Method 1: oracle JISST-PCA
        [u1_est, V1_est, W1_est, ~] = iMsst_sub(X{j}, Y{j}, u0, r, r, lambda, tol, max_iter);

        err_oracle{1, 1}(i, j) = abs(sin_d(u1_est{2}, u1)); % est err of u1
        err_oracle{1, 2}(i, j) = abs(sin_d(u1_est{3}, u2)); % est err of u2
        err_oracle{2, 1}(i, j) = abs(sin_d(V1_est{2}, V1)); % est err of V1
        err_oracle{2, 2}(i, j) = abs(sin_d(V1_est{3}, V2)); % est err of V2
        err_oracle{3, 1}(i, j) = abs(sin_d(W1_est{2}, W1)); % est err of W1
        err_oracle{3, 2}(i, j) = abs(sin_d(W1_est{3}, W2)); % est err of W2

        % Method 2: data-driven JisstPCA
        bic_r = bic_est(X{j}, Y{j}, p, q, 2, u0, lam, tol, max_iter); % estimated rank
        bic_r = bic_r';
        [u2_est, V2_est, W2_est, ~] = iMsst_sub(X{j}, Y{j}, u0, bic_r, bic_r, lambda, tol, max_iter);

        err_bic{1, 1}(i, j) = abs(sin_d(u2_est{2}, u1)); % est err of u1
        err_bic{1, 2}(i, j) = abs(sin_d(u2_est{3}, u2)); % est err of u2
        err_bic{2, 1}(i, j) = abs(sin_d(V2_est{2}, V1)); % est err of V1
        err_bic{2, 2}(i, j) = abs(sin_d(V2_est{3}, V2)); % est err of V2
        err_bic{3, 1}(i, j) = abs(sin_d(W2_est{2}, W1)); % est err of W1
        err_bic{3, 2}(i, j) = abs(sin_d(W2_est{3}, W2)); % est err of W2

        % Method 3: iHOSVD
        [u3_est, V3_est, W3_est, ~, ~] = iHOSVD(X{j}, Y{j}, r_sum, r_sum, 2);

        err_ihosvd{1, 1}(i, j) = abs(sin_d(u3_est(:, 1), u1));
        err_ihosvd{1, 2}(i, j) = abs(sin_d(u3_est(:, 2), u2));
        err_ihosvd{2, 1}(i, j) = abs(sin_d(V3_est(:, 1:cr(1)), V1));
        err_ihosvd{2, 2}(i, j) = abs(sin_d(V3_est(:, (cr(1)+1):cr(2)), V2));
        err_ihosvd{3, 1}(i, j) = abs(sin_d(W3_est(:, 1:cr(1)), W1));
        err_ihosvd{3, 2}(i, j) = abs(sin_d(W3_est(:, (cr(1)+1):cr(2)), W2));

        % Method 4: iHOOI
        [u4_est, V4_est, W4_est, ~, ~] = iHOOI(X{j}, Y{j}, r_sum, r_sum, max_iter, 2);

        err_ihooi{1, 1}(i, j) = abs(sin_d(u4_est(:, 1), u1));
        err_ihooi{1, 2}(i, j) = abs(sin_d(u4_est(:, 2), u2));
        err_ihooi{2, 1}(i, j) = abs(sin_d(V4_est(:, 1:cr(1)), V1));
        err_ihooi{2, 2}(i, j) = abs(sin_d(V4_est(:, (cr(1)+1):cr(2)), V2));
        err_ihooi{3, 1}(i, j) = abs(sin_d(W4_est(:, 1:cr(1)), W1));
        err_ihooi{3, 2}(i, j) = abs(sin_d(W4_est(:, (cr(1)+1):cr(2)), W2));
    end
end

% export the data to local file
writematrix(err_oracle{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_oracle_u1.csv")
writematrix(err_oracle{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_oracle_u2.csv")
writematrix(err_oracle{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_oracle_V1.csv")
writematrix(err_oracle{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_oracle_V2.csv")
writematrix(err_oracle{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_oracle_W1.csv")
writematrix(err_oracle{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_oracle_W2.csv")

writematrix(err_bic{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_bic_u1.csv")
writematrix(err_bic{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_bic_u2.csv")
writematrix(err_bic{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_bic_V1.csv")
writematrix(err_bic{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_bic_V2.csv")
writematrix(err_bic{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_bic_W1.csv")
writematrix(err_bic{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_bic_W2.csv")

writematrix(err_ihosvd{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihosvd_u1.csv")
writematrix(err_ihosvd{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihosvd_u2.csv")
writematrix(err_ihosvd{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihosvd_V1.csv")
writematrix(err_ihosvd{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihosvd_V2.csv")
writematrix(err_ihosvd{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihosvd_W1.csv")
writematrix(err_ihosvd{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihosvd_W2.csv")

writematrix(err_ihooi{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihooi_u1.csv")
writematrix(err_ihooi{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihooi_u2.csv")
writematrix(err_ihooi{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihooi_V1.csv")
writematrix(err_ihooi{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihooi_V2.csv")
writematrix(err_ihooi{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihooi_W1.csv")
writematrix(err_ihooi{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case1/err_ihooi_W2.csv")


%% scenario 2: p = q = 150, N = 50
p = 150;
q = 150;
N = 50;
r = [3 2];
dx = [2*p p];
dy = [2*q q];
tol = 0.0001;
max_iter = 100;

rep = 10;
eta = 0.4:0.4:2; % SNR level
sz = size(eta, 2);

err_oracle = cell(3, 2);
err_bic = cell(3, 2);
err_ihosvd = cell(3, 2);
err_ihooi = cell(3, 2);
for i = 1:3
    for j = 1:2
        err_oracle{i, j} = zeros(rep, sz);
        err_bic{i, j} = zeros(rep, sz);
        err_ihosvd{i, j} = zeros(rep, sz);
        err_ihooi{i, j} = zeros(rep, sz);
    end
end

rng("default")
V1 = orth(normrnd(0, 1, [p, r(1)]));
V2 = orth(normrnd(0, 1, [p, r(2)]));
W1 = orth(normrnd(0, 1, [q, r(1)]));
W2 = orth(normrnd(0, 1, [q, r(2)]));
u = normrnd(0, 1, [N, 3]);
u1 = u(:, 1)/norm(u(:, 1));
u2 = u(:, 2)/norm(u(:, 2));
X1 = squeeze(ttt(tensor(V1*V1'), tensor(u1)));
X2 = squeeze(ttt(tensor(V2*V2'), tensor(u2)));
Y1 = squeeze(ttt(tensor(W1*W1'), tensor(u1)));
Y2 = squeeze(ttt(tensor(W2*W2'), tensor(u2)));
X_t = dx(1)*X1 + dx(2)*X2;
Y_t = dy(1)*Y1 + dy(2)*Y2;

for i = 1:rep
    
    % set noise and perturbed tensor
    epi_X = zeros(size(X_t));
    epi_Y = zeros(size(Y_t));
    
    rng(i)
    for k = 1:N
        epi_X(:, :, k) = wigner(p, 1);
        epi_Y(:, :, k) = wigner(q, 1);
    end
    
    X = cell(sz, 1); 
    Y = cell(sz, 1);
    for j = 1:sz
        X{j} = tensor(eta(j)*(norm(tensor(epi_X))/norm(X_t))*X_t + epi_X);
        Y{j} = tensor(eta(j)*(norm(tensor(epi_Y))/norm(Y_t))*Y_t + epi_Y);
    end

    for j = 1:sz

        % initialization
        M = [double(tenmat(X{j}, 3, [1, 2])), double(tenmat(Y{j}, 3, [1, 2]))];
        [a, b, ~] = svd(M);
        [~, ind] = sort(diag(b));
        a = a(:, ind);
        a = fliplr(a);
        hat_u = a(:, 1);
        u0 = hat_u; % spectral initialization

        % scaler lambda
        lam = norm(X{j})/(norm(X{j}) + norm(Y{j}));
        lambda = [lam lam];

        % sum/cumsum of r
        r_sum = sum(r);
        cr = cumsum(r);

        % Method 1: oracle JISST-PCA
        [u1_est, V1_est, W1_est, ~] = iMsst_sub(X{j}, Y{j}, u0, r, r, lambda, tol, max_iter);

        err_oracle{1, 1}(i, j) = abs(sin_d(u1_est{2}, u1)); % est err of u1
        err_oracle{1, 2}(i, j) = abs(sin_d(u1_est{3}, u2)); % est err of u2
        err_oracle{2, 1}(i, j) = abs(sin_d(V1_est{2}, V1)); % est err of V1
        err_oracle{2, 2}(i, j) = abs(sin_d(V1_est{3}, V2)); % est err of V2
        err_oracle{3, 1}(i, j) = abs(sin_d(W1_est{2}, W1)); % est err of W1
        err_oracle{3, 2}(i, j) = abs(sin_d(W1_est{3}, W2)); % est err of W2

        % Method 2: data-driven JisstPCA
        bic_r = bic_est(X{j}, Y{j}, p, q, 2, u0, lam, tol, max_iter); % estimated rank
        bic_r = bic_r';
        [u2_est, V2_est, W2_est, ~] = iMsst_sub(X{j}, Y{j}, u0, bic_r, bic_r, lambda, tol, max_iter);

        err_bic{1, 1}(i, j) = abs(sin_d(u2_est{2}, u1)); % est err of u1
        err_bic{1, 2}(i, j) = abs(sin_d(u2_est{3}, u2)); % est err of u2
        err_bic{2, 1}(i, j) = abs(sin_d(V2_est{2}, V1)); % est err of V1
        err_bic{2, 2}(i, j) = abs(sin_d(V2_est{3}, V2)); % est err of V2
        err_bic{3, 1}(i, j) = abs(sin_d(W2_est{2}, W1)); % est err of W1
        err_bic{3, 2}(i, j) = abs(sin_d(W2_est{3}, W2)); % est err of W2

        % Method 3: iHOSVD
        [u3_est, V3_est, W3_est, ~, ~] = iHOSVD(X{j}, Y{j}, r_sum, r_sum, 2);

        err_ihosvd{1, 1}(i, j) = abs(sin_d(u3_est(:, 1), u1));
        err_ihosvd{1, 2}(i, j) = abs(sin_d(u3_est(:, 2), u2));
        err_ihosvd{2, 1}(i, j) = abs(sin_d(V3_est(:, 1:cr(1)), V1));
        err_ihosvd{2, 2}(i, j) = abs(sin_d(V3_est(:, (cr(1)+1):cr(2)), V2));
        err_ihosvd{3, 1}(i, j) = abs(sin_d(W3_est(:, 1:cr(1)), W1));
        err_ihosvd{3, 2}(i, j) = abs(sin_d(W3_est(:, (cr(1)+1):cr(2)), W2));

        % Method 4: iHOOI
        [u4_est, V4_est, W4_est, ~, ~] = iHOOI(X{j}, Y{j}, r_sum, r_sum, max_iter, 2);

        err_ihooi{1, 1}(i, j) = abs(sin_d(u4_est(:, 1), u1));
        err_ihooi{1, 2}(i, j) = abs(sin_d(u4_est(:, 2), u2));
        err_ihooi{2, 1}(i, j) = abs(sin_d(V4_est(:, 1:cr(1)), V1));
        err_ihooi{2, 2}(i, j) = abs(sin_d(V4_est(:, (cr(1)+1):cr(2)), V2));
        err_ihooi{3, 1}(i, j) = abs(sin_d(W4_est(:, 1:cr(1)), W1));
        err_ihooi{3, 2}(i, j) = abs(sin_d(W4_est(:, (cr(1)+1):cr(2)), W2));
    end
end

% export the data to local file
writematrix(err_oracle{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_oracle_u1.csv")
writematrix(err_oracle{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_oracle_u2.csv")
writematrix(err_oracle{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_oracle_V1.csv")
writematrix(err_oracle{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_oracle_V2.csv")
writematrix(err_oracle{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_oracle_W1.csv")
writematrix(err_oracle{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_oracle_W2.csv")

writematrix(err_bic{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_bic_u1.csv")
writematrix(err_bic{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_bic_u2.csv")
writematrix(err_bic{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_bic_V1.csv")
writematrix(err_bic{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_bic_V2.csv")
writematrix(err_bic{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_bic_W1.csv")
writematrix(err_bic{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_bic_W2.csv")

writematrix(err_ihosvd{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihosvd_u1.csv")
writematrix(err_ihosvd{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihosvd_u2.csv")
writematrix(err_ihosvd{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihosvd_V1.csv")
writematrix(err_ihosvd{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihosvd_V2.csv")
writematrix(err_ihosvd{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihosvd_W1.csv")
writematrix(err_ihosvd{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihosvd_W2.csv")

writematrix(err_ihooi{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihooi_u1.csv")
writematrix(err_ihooi{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihooi_u2.csv")
writematrix(err_ihooi{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihooi_V1.csv")
writematrix(err_ihooi{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihooi_V2.csv")
writematrix(err_ihooi{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihooi_W1.csv")
writematrix(err_ihooi{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case2/err_ihooi_W2.csv")

%% scenario 3: p = 150, q = 50, N = 200
p = 150;
q = 50;
N = 200;
r = [3 2];
dx = [2*p p];
dy = [2*q q];
tol = 0.0001;
max_iter = 100;

rep = 10;
eta = 0.4:0.4:2; % SNR level
sz = size(eta, 2);

err_oracle = cell(3, 2);
err_bic = cell(3, 2);
err_ihosvd = cell(3, 2);
err_ihooi = cell(3, 2);
for i = 1:3
    for j = 1:2
        err_oracle{i, j} = zeros(rep, sz);
        err_bic{i, j} = zeros(rep, sz);
        err_ihosvd{i, j} = zeros(rep, sz);
        err_ihooi{i, j} = zeros(rep, sz);
    end
end

rng("default")
V1 = orth(normrnd(0, 1, [p, r(1)]));
V2 = orth(normrnd(0, 1, [p, r(2)]));
W1 = orth(normrnd(0, 1, [q, r(1)]));
W2 = orth(normrnd(0, 1, [q, r(2)]));
u = normrnd(0, 1, [N, 3]);
u1 = u(:, 1)/norm(u(:, 1));
u2 = u(:, 2)/norm(u(:, 2));
X1 = squeeze(ttt(tensor(V1*V1'), tensor(u1)));
X2 = squeeze(ttt(tensor(V2*V2'), tensor(u2)));
Y1 = squeeze(ttt(tensor(W1*W1'), tensor(u1)));
Y2 = squeeze(ttt(tensor(W2*W2'), tensor(u2)));
X_t = dx(1)*X1 + dx(2)*X2;
Y_t = dy(1)*Y1 + dy(2)*Y2;

for i = 1:rep
    
    % set noise and perturbed tensor
    epi_X = zeros(size(X_t));
    epi_Y = zeros(size(Y_t));
    
    rng(i)
    for k = 1:N
        epi_X(:, :, k) = wigner(p, 1);
        epi_Y(:, :, k) = wigner(q, 1);
    end
    
    X = cell(sz, 1); 
    Y = cell(sz, 1);
    for j = 1:sz
        X{j} = tensor(eta(j)*(norm(tensor(epi_X))/norm(X_t))*X_t + epi_X);
        Y{j} = tensor(eta(j)*(norm(tensor(epi_Y))/norm(Y_t))*Y_t + epi_Y);
    end

    for j = 1:sz

        % initialization
        M = [double(tenmat(X{j}, 3, [1, 2])), double(tenmat(Y{j}, 3, [1, 2]))];
        [a, b, ~] = svd(M);
        [~, ind] = sort(diag(b));
        a = a(:, ind);
        a = fliplr(a);
        hat_u = a(:, 1);
        u0 = hat_u; % spectral initialization

        % scaler lambda
        lam = norm(X{j})/(norm(X{j}) + norm(Y{j}));
        lambda = [lam lam];

        % sum/cumsum of r
        r_sum = sum(r);
        cr = cumsum(r);

        % Method 1: oracle JISST-PCA
        [u1_est, V1_est, W1_est, ~] = iMsst_sub(X{j}, Y{j}, u0, r, r, lambda, tol, max_iter);

        err_oracle{1, 1}(i, j) = abs(sin_d(u1_est{2}, u1)); % est err of u1
        err_oracle{1, 2}(i, j) = abs(sin_d(u1_est{3}, u2)); % est err of u2
        err_oracle{2, 1}(i, j) = abs(sin_d(V1_est{2}, V1)); % est err of V1
        err_oracle{2, 2}(i, j) = abs(sin_d(V1_est{3}, V2)); % est err of V2
        err_oracle{3, 1}(i, j) = abs(sin_d(W1_est{2}, W1)); % est err of W1
        err_oracle{3, 2}(i, j) = abs(sin_d(W1_est{3}, W2)); % est err of W2

        % Method 2: data-driven JisstPCA
        bic_r = bic_est(X{j}, Y{j}, p, q, 2, u0, lam, tol, max_iter); % estimated rank
        bic_r = bic_r';
        [u2_est, V2_est, W2_est, ~] = iMsst_sub(X{j}, Y{j}, u0, bic_r, bic_r, lambda, tol, max_iter);

        err_bic{1, 1}(i, j) = abs(sin_d(u2_est{2}, u1)); % est err of u1
        err_bic{1, 2}(i, j) = abs(sin_d(u2_est{3}, u2)); % est err of u2
        err_bic{2, 1}(i, j) = abs(sin_d(V2_est{2}, V1)); % est err of V1
        err_bic{2, 2}(i, j) = abs(sin_d(V2_est{3}, V2)); % est err of V2
        err_bic{3, 1}(i, j) = abs(sin_d(W2_est{2}, W1)); % est err of W1
        err_bic{3, 2}(i, j) = abs(sin_d(W2_est{3}, W2)); % est err of W2

        % Method 3: iHOSVD
        [u3_est, V3_est, W3_est, ~, ~] = iHOSVD(X{j}, Y{j}, r_sum, r_sum, 2);

        err_ihosvd{1, 1}(i, j) = abs(sin_d(u3_est(:, 1), u1));
        err_ihosvd{1, 2}(i, j) = abs(sin_d(u3_est(:, 2), u2));
        err_ihosvd{2, 1}(i, j) = abs(sin_d(V3_est(:, 1:cr(1)), V1));
        err_ihosvd{2, 2}(i, j) = abs(sin_d(V3_est(:, (cr(1)+1):cr(2)), V2));
        err_ihosvd{3, 1}(i, j) = abs(sin_d(W3_est(:, 1:cr(1)), W1));
        err_ihosvd{3, 2}(i, j) = abs(sin_d(W3_est(:, (cr(1)+1):cr(2)), W2));

        % Method 4: iHOOI
        [u4_est, V4_est, W4_est, ~, ~] = iHOOI(X{j}, Y{j}, r_sum, r_sum, max_iter, 2);

        err_ihooi{1, 1}(i, j) = abs(sin_d(u4_est(:, 1), u1));
        err_ihooi{1, 2}(i, j) = abs(sin_d(u4_est(:, 2), u2));
        err_ihooi{2, 1}(i, j) = abs(sin_d(V4_est(:, 1:cr(1)), V1));
        err_ihooi{2, 2}(i, j) = abs(sin_d(V4_est(:, (cr(1)+1):cr(2)), V2));
        err_ihooi{3, 1}(i, j) = abs(sin_d(W4_est(:, 1:cr(1)), W1));
        err_ihooi{3, 2}(i, j) = abs(sin_d(W4_est(:, (cr(1)+1):cr(2)), W2));
    end
end

% export the data to local file
writematrix(err_oracle{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_oracle_u1.csv")
writematrix(err_oracle{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_oracle_u2.csv")
writematrix(err_oracle{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_oracle_V1.csv")
writematrix(err_oracle{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_oracle_V2.csv")
writematrix(err_oracle{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_oracle_W1.csv")
writematrix(err_oracle{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_oracle_W2.csv")

writematrix(err_bic{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_bic_u1.csv")
writematrix(err_bic{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_bic_u2.csv")
writematrix(err_bic{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_bic_V1.csv")
writematrix(err_bic{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_bic_V2.csv")
writematrix(err_bic{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_bic_W1.csv")
writematrix(err_bic{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_bic_W2.csv")

writematrix(err_ihosvd{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihosvd_u1.csv")
writematrix(err_ihosvd{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihosvd_u2.csv")
writematrix(err_ihosvd{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihosvd_V1.csv")
writematrix(err_ihosvd{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihosvd_V2.csv")
writematrix(err_ihosvd{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihosvd_W1.csv")
writematrix(err_ihosvd{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihosvd_W2.csv")

writematrix(err_ihooi{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihooi_u1.csv")
writematrix(err_ihooi{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihooi_u2.csv")
writematrix(err_ihooi{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihooi_V1.csv")
writematrix(err_ihooi{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihooi_V2.csv")
writematrix(err_ihooi{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihooi_W1.csv")
writematrix(err_ihooi{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case3/err_ihooi_W2.csv")

%% scenario 4: p = 150, q = 50, N = 50
p = 150;
q = 50;
N = 50;
r = [3 2];
dx = [2*p p];
dy = [2*q q];
tol = 0.0001;
max_iter = 100;

rep = 10;
eta = 0.4:0.4:2; % SNR level
sz = size(eta, 2);

err_oracle = cell(3, 2);
err_bic = cell(3, 2);
err_ihosvd = cell(3, 2);
err_ihooi = cell(3, 2);
for i = 1:3
    for j = 1:2
        err_oracle{i, j} = zeros(rep, sz);
        err_bic{i, j} = zeros(rep, sz);
        err_ihosvd{i, j} = zeros(rep, sz);
        err_ihooi{i, j} = zeros(rep, sz);
    end
end

rng("default")
V1 = orth(normrnd(0, 1, [p, r(1)]));
V2 = orth(normrnd(0, 1, [p, r(2)]));
W1 = orth(normrnd(0, 1, [q, r(1)]));
W2 = orth(normrnd(0, 1, [q, r(2)]));
u = normrnd(0, 1, [N, 3]);
u1 = u(:, 1)/norm(u(:, 1));
u2 = u(:, 2)/norm(u(:, 2));
X1 = squeeze(ttt(tensor(V1*V1'), tensor(u1)));
X2 = squeeze(ttt(tensor(V2*V2'), tensor(u2)));
Y1 = squeeze(ttt(tensor(W1*W1'), tensor(u1)));
Y2 = squeeze(ttt(tensor(W2*W2'), tensor(u2)));
X_t = dx(1)*X1 + dx(2)*X2;
Y_t = dy(1)*Y1 + dy(2)*Y2;

for i = 1:rep
    
    % set noise and perturbed tensor
    epi_X = zeros(size(X_t));
    epi_Y = zeros(size(Y_t));
    
    rng(i)
    for k = 1:N
        epi_X(:, :, k) = wigner(p, 1);
        epi_Y(:, :, k) = wigner(q, 1);
    end
    
    X = cell(sz, 1); 
    Y = cell(sz, 1);
    for j = 1:sz
        X{j} = tensor(eta(j)*(norm(tensor(epi_X))/norm(X_t))*X_t + epi_X);
        Y{j} = tensor(eta(j)*(norm(tensor(epi_Y))/norm(Y_t))*Y_t + epi_Y);
    end

    for j = 1:sz

        % initialization
        M = [double(tenmat(X{j}, 3, [1, 2])), double(tenmat(Y{j}, 3, [1, 2]))];
        [a, b, ~] = svd(M);
        [~, ind] = sort(diag(b));
        a = a(:, ind);
        a = fliplr(a);
        hat_u = a(:, 1);
        u0 = hat_u; % spectral initialization

        % scaler lambda
        lam = norm(X{j})/(norm(X{j}) + norm(Y{j}));
        lambda = [lam lam];

        % sum/cumsum of r
        r_sum = sum(r);
        cr = cumsum(r);

        % Method 1: oracle JISST-PCA
        [u1_est, V1_est, W1_est, ~] = iMsst_sub(X{j}, Y{j}, u0, r, r, lambda, tol, max_iter);

        err_oracle{1, 1}(i, j) = abs(sin_d(u1_est{2}, u1)); % est err of u1
        err_oracle{1, 2}(i, j) = abs(sin_d(u1_est{3}, u2)); % est err of u2
        err_oracle{2, 1}(i, j) = abs(sin_d(V1_est{2}, V1)); % est err of V1
        err_oracle{2, 2}(i, j) = abs(sin_d(V1_est{3}, V2)); % est err of V2
        err_oracle{3, 1}(i, j) = abs(sin_d(W1_est{2}, W1)); % est err of W1
        err_oracle{3, 2}(i, j) = abs(sin_d(W1_est{3}, W2)); % est err of W2

        % Method 2: data-driven JisstPCA
        bic_r = bic_est(X{j}, Y{j}, p, q, 2, u0, lam, tol, max_iter); % estimated rank
        bic_r = bic_r';
        [u2_est, V2_est, W2_est, ~] = iMsst_sub(X{j}, Y{j}, u0, bic_r, bic_r, lambda, tol, max_iter);

        err_bic{1, 1}(i, j) = abs(sin_d(u2_est{2}, u1)); % est err of u1
        err_bic{1, 2}(i, j) = abs(sin_d(u2_est{3}, u2)); % est err of u2
        err_bic{2, 1}(i, j) = abs(sin_d(V2_est{2}, V1)); % est err of V1
        err_bic{2, 2}(i, j) = abs(sin_d(V2_est{3}, V2)); % est err of V2
        err_bic{3, 1}(i, j) = abs(sin_d(W2_est{2}, W1)); % est err of W1
        err_bic{3, 2}(i, j) = abs(sin_d(W2_est{3}, W2)); % est err of W2

        % Method 3: iHOSVD
        [u3_est, V3_est, W3_est, ~, ~] = iHOSVD(X{j}, Y{j}, r_sum, r_sum, 2);

        err_ihosvd{1, 1}(i, j) = abs(sin_d(u3_est(:, 1), u1));
        err_ihosvd{1, 2}(i, j) = abs(sin_d(u3_est(:, 2), u2));
        err_ihosvd{2, 1}(i, j) = abs(sin_d(V3_est(:, 1:cr(1)), V1));
        err_ihosvd{2, 2}(i, j) = abs(sin_d(V3_est(:, (cr(1)+1):cr(2)), V2));
        err_ihosvd{3, 1}(i, j) = abs(sin_d(W3_est(:, 1:cr(1)), W1));
        err_ihosvd{3, 2}(i, j) = abs(sin_d(W3_est(:, (cr(1)+1):cr(2)), W2));

        % Method 4: iHOOI
        [u4_est, V4_est, W4_est, ~, ~] = iHOOI(X{j}, Y{j}, r_sum, r_sum, max_iter, 2);

        err_ihooi{1, 1}(i, j) = abs(sin_d(u4_est(:, 1), u1));
        err_ihooi{1, 2}(i, j) = abs(sin_d(u4_est(:, 2), u2));
        err_ihooi{2, 1}(i, j) = abs(sin_d(V4_est(:, 1:cr(1)), V1));
        err_ihooi{2, 2}(i, j) = abs(sin_d(V4_est(:, (cr(1)+1):cr(2)), V2));
        err_ihooi{3, 1}(i, j) = abs(sin_d(W4_est(:, 1:cr(1)), W1));
        err_ihooi{3, 2}(i, j) = abs(sin_d(W4_est(:, (cr(1)+1):cr(2)), W2));
    end
end

% export the data to local file
writematrix(err_oracle{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_oracle_u1.csv")
writematrix(err_oracle{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_oracle_u2.csv")
writematrix(err_oracle{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_oracle_V1.csv")
writematrix(err_oracle{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_oracle_V2.csv")
writematrix(err_oracle{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_oracle_W1.csv")
writematrix(err_oracle{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_oracle_W2.csv")

writematrix(err_bic{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_bic_u1.csv")
writematrix(err_bic{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_bic_u2.csv")
writematrix(err_bic{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_bic_V1.csv")
writematrix(err_bic{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_bic_V2.csv")
writematrix(err_bic{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_bic_W1.csv")
writematrix(err_bic{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_bic_W2.csv")

writematrix(err_ihosvd{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihosvd_u1.csv")
writematrix(err_ihosvd{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihosvd_u2.csv")
writematrix(err_ihosvd{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihosvd_V1.csv")
writematrix(err_ihosvd{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihosvd_V2.csv")
writematrix(err_ihosvd{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihosvd_W1.csv")
writematrix(err_ihosvd{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihosvd_W2.csv")

writematrix(err_ihooi{1, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihooi_u1.csv")
writematrix(err_ihooi{1, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihooi_u2.csv")
writematrix(err_ihooi{2, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihooi_V1.csv")
writematrix(err_ihooi{2, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihooi_V2.csv")
writematrix(err_ihooi{3, 1}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihooi_W1.csv")
writematrix(err_ihooi{3, 2}, "/Users/jl259/Desktop/section_4/unstructure/case4/err_ihooi_W2.csv")



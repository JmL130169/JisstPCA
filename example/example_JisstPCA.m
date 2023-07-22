% An illustrative example of JisstPCA

%% data setup
p = 30; q = 20; N = 80; % dimension of tensors
K = 2;
rx = [3 2]; ry = rx;
dx = [2*p p]; % signal of X for each layer
dy = [2*q q]; % signal of Y for each layer

rng(2023)
V = cell(K, 1); W = cell(K, 1); u = cell(K, 1);
X = zeros(p, p, N); Y = zeros(q, q, N);
for k = 1 : K
    V{k} = orth(normrnd(0, 1, [p, rx(k)])); % kth layer true factor V1
    W{k} = orth(normrnd(0, 1, [q, ry(k)])); % kth layer true factor W1
    u{k} = orth(normrnd(0, 1, [N, 1]));
    X = X + dx(k) * squeeze(ttt(tensor(V{k}*V{k}'), tensor(u{k}))); % noiseless X
    Y = Y + dy(k) * squeeze(ttt(tensor(W{k}*W{k}'), tensor(u{k}))); % noiseless Y
end

noise_X = zeros(size(X));
noise_Y = zeros(size(Y));
for k = 1:N
    noise_X(:, :, k) = wigner(p, 1); 
    noise_Y(:, :, k) = wigner(q, 1); % generate noise from Gaussian ensemble, by each slice
end

X_obs = tensor(X + noise_X); % observation tensor of X, with noise
Y_obs = tensor(Y + noise_Y); % observation tensor of X, with noise

% optional inputs
tol = 0.0001; % tolerance level
max_iter = 100; % maximum iteration number
ratio = norm(X_obs)/(norm(X_obs) + norm(Y_obs)); % relative weight of X and Y in estimation
lambda = ratio*ones(1, K); % relative weight of X and Y
u0 = init(X_obs, Y_obs, lambda(1)); % initialization
deflation = 0;

%% case 1: only X and Y are available
K = 2;
[u_est1, V_est1, W_est1, d_est1] = JisstPCA(X_obs, Y_obs, K);


%% case 2: rx and ry are also available
[u_est2, V_est2, W_est2, d_est2] = JisstPCA(X_obs, Y_obs, K, 'rx', rx, 'ry', ry);


%% case 3: all the inputs are available
[u_est3, V_est3, W_est3, d_est3] = JisstPCA(X_obs, Y_obs, K, 'rx', rx, 'ry', ry, ...
    'lambda', lambda, 'u0', u0, 'tol', tol, 'max_iter', max_iter, ...
    'deflation', deflation);
find_dist(u_est3, V_est3, W_est3, u, V, W)






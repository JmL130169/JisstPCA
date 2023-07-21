% An illustrative example of JisstPCA

%% data setup
p = 30; q = 20; N = 80; % dimension of tensors
rx = [3 2]; ry = rx;
dx = [2*p p]; % signal of X for each layer
dy = [2*q q]; % signal of Y for each layer

rng(2023)
V1 = orth(normrnd(0, 1, [p, rx(1)])); % first layer true factor V1
W1 = orth(normrnd(0, 1, [q, ry(1)])); % first layer true factor W1
V2 = orth(normrnd(0, 1, [p, rx(2)])); % first layer true factor V2
W2 = orth(normrnd(0, 1, [q, ry(2)])); % first layer true factor W2
u = normrnd(0, 1, [N, 2]);
u1 = u(:, 1)/norm(u(:, 1)); % first joint factor u1
u2 = u(:, 2)/norm(u(:, 2)); % second joint factor u2

X = dx(1)*squeeze(ttt(tensor(V1*V1'), tensor(u1))) + dx(2)*squeeze(ttt(tensor(V2*V2'), tensor(u2))); % noiseless X
Y = dy(1)*squeeze(ttt(tensor(W1*W1'), tensor(u1))) + dy(2)*squeeze(ttt(tensor(W2*W2'), tensor(u2))); % noiseless Y

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
lambda = ratio*ones(1, 2); % relative weight of X and Y
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
norm(u_est3{1} * u_est3{1}' - u1 * u1')
norm(u_est3{2} * u_est3{2}' - u2 * u2')
norm(V_est3{1} * V_est3{1}' - V1 * V1')
norm(V_est3{2} * V_est3{2}' - V2 * V2')
norm(W_est3{1} * W_est3{1}' - W1 * W1')
norm(W_est3{2} * W_est3{2}' - W2 * W2')






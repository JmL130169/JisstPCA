% This is an toy example of single-factor/multi-factor JisstPCA
% The signals are represented in scaler dx and dy, so we apply Jisst_single and
% Jisst_multi

%% single-factor JisstPCA

% data setup

p = 50; 
q = 50; 
N = 100; 
rx = 3;
ry = 2;
dx = 5*p; % signal of X
dy = 4*q; % signal of Y

rng(2023)
V = orth(normrnd(0, 1, [p, rx])); % true factor V, generated from iid N(0, 1)
W = orth(normrnd(0, 1, [q, ry])); % true factor W
u = normrnd(0, 1, [N, 1]);
u = u/norm(u); % normalize joint factor u

X = dx*squeeze(ttt(tensor(V*V'), tensor(u))); % noiseless X
Y = dy*squeeze(ttt(tensor(W*W'), tensor(u))); % noiseless Y

noise_X = zeros(size(X));
noise_Y = zeros(size(Y));
for k = 1:N
    noise_X(:, :, k) = wigner(p, 1); 
    noise_Y(:, :, k) = wigner(q, 1); % generate noise from Gaussian ensemble, by each slice
end

X_obs = tensor(X + noise_X); % observation tensor of X, with noise
Y_obs = tensor(Y + noise_Y); % observation tensor of X, with noise


% estimate tensor factors V, W, u with noised data, using Jisst_single

u0 = init(X_obs, Y_obs); % initialization
tol = 0.0001; % tolerance level
max_iter = 100; % maximum iteration number
lambda = norm(X_obs)/(norm(X_obs) + norm(Y_obs)); % relative weight of X and Y in estimation

[hat_u, hat_V, hat_W, d_x, d_y, ~, ~] = Jisst_single(X_obs, Y_obs, u0, rx, ry, lambda, tol, max_iter);
% by this step, the output is estimation of each tensor factor

%% multi-factor JisstPCA, K = 2

% data setup

p = 50; 
q = 50; 
N = 100; 
rx = [4 3];
ry = [3 2];
dx = [10*p 5*p]; % signal of X for each layer
dy = [8*q 4*q]; % signal of Y for each layer

rng(2023)
V1 = orth(normrnd(0, 1, [p, rx(1)])); % first layer true factor V1
W1 = orth(normrnd(0, 1, [q, ry(1)])); % first layer true factor W1
V2 = orth(normrnd(0, 1, [p, rx(2)])); % first layer true factor V2
W2 = orth(normrnd(0, 1, [q, ry(2)])); % first layer true factor W2
u = orth(normrnd(0, 1, [N, 2]));
u1 = u(:, 1); % first joint factor u1
u2 = u(:, 2); % second joint factor u2

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


% estimate tensor factors Vk, Wk, uk with noised data, using Jisst_multi

u0 = init(X_obs, Y_obs); % initialization
tol = 0.0001; % tolerance level
max_iter = 100; % maximum iteration number
ratio = norm(X_obs)/(norm(X_obs) + norm(Y_obs)); % relative weight of X and Y in estimation
lambda = ratio*ones([1, 2]); % relative weight of X and Y

[u_est, V_est, W_est, d_est] = Jisst_multi(X_obs, Y_obs, u0, rx, ry, lambda, tol, max_iter);
% by this step, we get the estimation of each tensor factor
% Note: as mentioned in the function, the estimation of uk, Vk, Wk are
% u_est{k+1}, V_est{k+1} and W_est{k+1}, respectively. For 2-(K+1) matrix
% d_est, we have d_est(1, 2:K+1) as the estimation of dx and d_est(2,
% 2:K+1) as the estimation of dy










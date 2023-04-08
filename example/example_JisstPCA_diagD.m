% This is an toy example of single-factor/multi-factor JisstPCA
% The signals are represented in diagonal matrices Dx and Dy, so we apply dJisst_single and
% dJisst_multi

%% single-factor JisstPCA, diagonal signal D

% data setup

p = 50; 
q = 50; 
N = 100; 
rx = 3;
ry = 2;
Dx = diag([10*p 8*p 6*p]); % signal of X, diagonal matrix Dx
Dy = diag([8*q 6*q]); % signal of Y, diagonal matrix Dy

rng(2023)
V = orth(normrnd(0, 1, [p, rx])); % true factor V, generated from iid N(0, 1)
W = orth(normrnd(0, 1, [q, ry])); % true factor W
u = normrnd(0, 1, [N, 1]);
u = u/norm(u); % normalize joint factor u

X = squeeze(ttt(tensor(V*Dx*V'), tensor(u))); % noiseless X
Y = squeeze(ttt(tensor(W*Dy*W'), tensor(u))); % noiseless Y

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

[hat_u, hat_V, hat_W, D_x, D_y, ~, ~] = dJisst_single(X_obs, Y_obs, u0, rx, ry, lambda, tol, max_iter);
% output is estimation of each tensor factor

%% multi-factor JisstPCA, K = 2, diagonal signal D

% data setup

p = 50; 
q = 50; 
N = 100; 
rx = [4 3];
ry = [3 2];
Dx1 = diag([10*p 8*p 6*p 4*p]); % signal of X, first layer
Dy1 = diag([8*q 6*q 4*q]); % signal of Y, first layer
Dx2 = diag([3*p 2*p p]); % signal of X, second layer
Dy2 = diag([3*q 2*q]); % signal of Y, second layer

rng(2023)
V1 = orth(normrnd(0, 1, [p, rx(1)])); % first layer true factor V1
W1 = orth(normrnd(0, 1, [q, ry(1)])); % first layer true factor W1
V2 = orth(normrnd(0, 1, [p, rx(2)])); % first layer true factor V2
W2 = orth(normrnd(0, 1, [q, ry(2)])); % first layer true factor W2
u = orth(normrnd(0, 1, [N, 2]));
u1 = u(:, 1); % first joint factor u1
u2 = u(:, 2); % second joint factor u2

X = squeeze(ttt(tensor(V1*Dx1*V1'), tensor(u1))) + squeeze(ttt(tensor(V2*Dx2*V2'), tensor(u2))); % noiseless X
Y = squeeze(ttt(tensor(W1*Dy1*W1'), tensor(u1))) + squeeze(ttt(tensor(W2*Dy2*W2'), tensor(u2))); % noiseless Y

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

[u_est, V_est, W_est, Dx_est, Dy_est] = dJisst_multi(X_obs, Y_obs, u0, rx, ry, lambda, tol, max_iter);
% by this step, we get the estimation of each tensor factor
% Note: as mentioned in the function, the estimation of uk, Vk, Wk are
% u_est{k+1}, V_est{k+1} and W_est{k+1}, respectively. And Dxk, Dyk are
% estimated by Dx_est{k+1}, Dy_est{k+1}, respectively

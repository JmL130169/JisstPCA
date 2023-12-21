%% theory validation for single-factor case
% Software dependency: tensor toolbox (https://www.tensortoolbox.org)
addpath('JisstPCA/functions')
filename = "JisstPCA/simulations/TheoryValidation/JisstPCA_res.csv";

%% data construction for validating theory with scalar d
%model setting
p_vec = [60 90 120];
SNR_scale_vec = [1.5 3 6 12 24];
rx = 3;
ry = 2;

%default parameters
lam = 0.5;

num_p = length(p_vec);
num_SNR = length(SNR_scale_vec);
nrep = 20; % repeats
max_iter = 10;

ii = 0;
for l = 1 : num_p
    p = p_vec(l);
    q = p;
    N = 2 * p;

    %generate data
    rng(1234)
    u = normrnd(0, 1, [N, 1]);
    u = u/norm(u);
    V = orth(normrnd(0, 1, [p, rx])); % true factor V
    W = orth(normrnd(0, 1, [q, ry])); % true factor W
    for j = 1 : num_SNR
        SNR_scale = SNR_scale_vec(j);
        dx = SNR_scale * (sqrt(p) + sqrt(N));
        dy = 1.2 * dx;
        X = dx * squeeze(ttt(tensor(V*V'), tensor(u)));
        Y = dy * squeeze(ttt(tensor(W*W'), tensor(u)));
        for k = 1 : nrep
            ii = ii + 1;
            rng(k)
            epi_X = wigner_sst(p, N);
            epi_Y = wigner_sst(q, N);
            X_obs = tensor(double(X) + epi_X);
            Y_obs = tensor(double(Y) + epi_Y);

            %estimation
            M = [double(lam * tenmat(X_obs, 3, [1, 2])), double((1 - lam) * tenmat(Y_obs, 3, [1, 2]))];
            [u0, ~, ~] = svds(M, 1);
            [hat_u, hat_V, hat_W, loss] = Jisst_single_iter(X_obs, Y_obs, u0, rx, ry, lam, max_iter, true);
            for kk = 0 : max_iter
                if kk == 0
                    factor = "u"; dim = p; SNR = SNR_scale;
                    seed = k; iteration = kk;
                    stat_err = abs(sin_do(hat_u{kk + 1}, u));
                    optim_loss = NaN;
                else
                    factor = ["u"; "V"; "W"]; dim = p * ones(3, 1); SNR = SNR_scale * ones(3, 1);
                    seed = k * ones(3, 1); iteration = kk * ones(3, 1);
                    stat_err = [abs(sin_do(hat_u{kk + 1}, u)); abs(sin_do(hat_V{kk}, V)); abs(sin_do(hat_W{kk}, W))];
                    optim_loss = loss(kk) * ones(3, 1);
                end
                results = table(factor, dim, SNR, seed, iteration, stat_err, optim_loss);
                if isfile(filename)
                    writetable(results, filename,...
                    'WriteMode','Append','WriteVariableNames',false);
                else
                    writetable(results, filename);
                end
            end
            sprintf("finished p = %d, SNR = %d, seed = %d", p, SNR_scale, k)
        end
    end
end


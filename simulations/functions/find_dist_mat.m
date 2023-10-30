function results = find_dist_mat(u_est, V_est, W_est, u, V, W)
    K = length(u);
    results = zeros(3, K);
    start_V = 1; start_W = 1;
    for k = 1 : K
        rx = size(V{k}, 2); ry = size(W{k}, 2);
        results(1, k) = dist(u_est(:, k), u{k});
        results(2, k) = dist(V_est(:, start_V : (start_V + rx - 1)), V{k});
        results(3, k) = dist(W_est(:, start_W : (start_W + ry - 1)), W{k});
        start_V = start_V + rx; start_W = start_W + ry;
    end
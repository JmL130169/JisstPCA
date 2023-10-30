function results = find_dist(u_est, V_est, W_est, u, V, W)
    K = length(u_est);
    results = zeros(3, K);
    for k = 1 : K
        results(1, k) = dist(u_est{k}, u{k});
        results(2, k) = dist(V_est{k}, V{k});
        results(3, k) = dist(W_est{k}, W{k});
    end
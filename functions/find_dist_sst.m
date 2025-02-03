function results = find_dist_sst(ux_est, uy_est, V_est, W_est, ux, uy, V, W)

    K = length(ux_est);
    results = zeros(4, K);
    for k = 1 : K
        results(1, k) = dist(ux_est{k}, ux{k});
        results(2, k) = dist(uy_est{k}, uy{k});
        results(3, k) = dist(V_est{k}, V{k});
        results(4, k) = dist(W_est{k}, W{k});
    end

end
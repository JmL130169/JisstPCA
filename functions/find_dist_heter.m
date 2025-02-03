function results = find_dist_heter(u_est, V_est, W_est, ux, uy, V, W)

    K = length(u_est);
    results = zeros(4, K);
    for k = 1 : K
        results(1, k) = dist(u_est{k}, ux{k});
        results(2, k) = dist(u_est{k}, uy{k});
        results(3, k) = dist(V_est{k}, V{k});
        results(4, k) = dist(W_est{k}, W{k});
    end

end
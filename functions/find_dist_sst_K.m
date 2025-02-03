function results = find_dist_sst_K(ux_est, uy_est, V_est, W_est, ux, uy, V, W)

    K = length(ux_est);
    Kx = length(V_est);
    Ky = length(W_est);
    results = zeros(4, K);

    for k = 1 : Kx
        results(1, k) = dist(ux_est{k}, ux{k});
        results(3, k) = dist(V_est{k}, V{k});
    end

    for k = 1 : Ky
        results(2, k) = dist(uy_est{k}, uy{k});
        results(4, k) = dist(W_est{k}, W{k});
    end

end
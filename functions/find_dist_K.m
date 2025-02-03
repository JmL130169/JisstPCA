function results = find_dist_K(u_est, V_est, W_est, u, V, W)

    Kx = length(V);
    Ky = length(W);
    K = max(Kx, Ky);
    results = zeros(3, K);

    for k = 1 : K
        results(1, k) = dist(u_est{k}, u{k});
    end

    for k = 1 : Kx
        results(2, k) = dist(V_est{k}, V{k});
    end

    for k = 1 : Ky
        results(3, k) = dist(W_est{k}, W{k});
    end

end
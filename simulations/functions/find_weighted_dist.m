function results = find_weighted_dist(u_est, V_est, Dx_est, W_est, Dy_est, u, Network_x, Network_y)
    K = length(u_est);
    results = zeros(3, K);
    for k = 1 : K
        results(1, k) = dist(u_est{k}, u{k});
        if(size(Dx_est, 1) == 1)
            results(2, k) = norm(V_est{k} * Dx_est(k) *  V_est{k}' - Network_x{k})/norm(Network_x{k});
            results(3, k) = norm(W_est{k} * Dy_est(k) * W_est{k}' - Network_y{k})/norm(Network_y{k});
        else
            results(2, k) = norm(V_est{k} * Dx_est{k} * V_est{k}' - Network_x{k})/norm(Network_x{k});
            results(3, k) = norm(W_est{k} * Dy_est{k} * W_est{k}' - Network_y{k})/norm(Network_y{k});
        end
    end
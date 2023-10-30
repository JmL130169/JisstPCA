function est_r = find_ranks(V_est, K)
    est_r = zeros(1, K); 
    for k = 1 : K
        est_r(k) = size(V_est{k}, 2);
    end
end
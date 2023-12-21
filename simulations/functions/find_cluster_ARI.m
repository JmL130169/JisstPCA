function [RI, ARI] = find_cluster_ARI(u_est, V_est, W_est, Dx_est, Dy_est, q, u, blocks, pathname)

    K = size(u_est, 1);
    RI = zeros(2 * K + 1, 1); ARI = zeros(2 * K + 1, 1);
    u_mat = u_est{1};
    for k = 2:K
        u_mat = cat(2, u_mat, u_est{k});
    end
    est_cluster = kmeans(u_mat, K);
    true_cluster = (u{1} == 0) + 1;
    [RI(1), ARI(1)] = randindex(est_cluster, true_cluster);

    for k = 1:K
        est_clusterX = kmeans(V_est{k} * Dx_est{k}, blocks(1, k));
        est_clusterY = kmeans(W_est{k} * Dy_est{k}, blocks(2, k));
        true_clusterX = table2array(readtable(strcat(pathname,sprintf("blocknumber_V%d_q%d.csv", k, q))));
        true_clusterY = table2array(readtable(strcat(pathname, sprintf("blocknumber_W%d_q%d.csv", k, q))));
        [RI(2 + 2 * (k - 1)), ARI(2 + 2 * (k - 1))] = randindex(est_clusterX, true_clusterX);
        [RI(3 + 2 * (k - 1)), ARI(3 + 2 * (k - 1))] = randindex(est_clusterY, true_clusterY);
    end
end
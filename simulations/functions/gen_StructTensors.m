function [X, Y, u, V, W] = gen_StructTensors(p, q, N, rx, ry, K, structure_list, Dx, Dy, sample_clusters)
    u = cell(K, 1); V = cell(K, 1); W = cell(K, 1);
    cluster_member = ceil(rand(N, 1) * sample_clusters);
    cluster_mean = rand(K, sample_clusters);% - 0.5;
    X = tensor(zeros(p, p, N)); Y = tensor(zeros(q, q, N)); 
    for k = 1 : K
        u{k} = orth(normrnd(0, 0.05, [N 1]) + cluster_mean(k, cluster_member)');
        block_V = sort(ceil(rand(p, 1) * rx(k)));
        block_W = sort(ceil(rand(q, 1) * ry(k)));
        if structure_list(k) == "block"
            V{k} = zeros(p, rx(k)); W{k} = zeros(q, ry(k));
            for j = 1 : rx(k)
                V{k}(block_V == j, j) = orth(normrnd(3, 1, [sum(block_V == j), 1]));
            end
            for j = 1 : ry(k)
                W{k}(block_W == j, j) = orth(normrnd(3, 1, [sum(block_W == j), 1]));
            end
        else if structure_list(k) == "star"
            for j = 1 : rx(k)
                graph = zeros(sum(block_V == j), sum(block_V == j));
                graph(:, 1) = 1; graph(1, :) = 1; graph(1, 1) = 0; laplacian = diag(sum(graph)) - graph;
                [V{k}(block_V == j, j), ~, ~] = svds(laplacian, 1);
            end
            for j = 1 : ry(k)
                graph = zeros(sum(block_W == j), sum(block_W == j)); 
                graph(:, 1) = 1; graph(1, :) = 1; graph(1, 1) = 0; laplacian = diag(sum(graph)) - graph;
                [W{k}(block_W == j, j), ~, ~] = svds(laplacian, 1);
            end
        else
            disp("Structure names have to be block or star!")
        end
        end
        X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(u{k})));
        Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(u{k})));
    end
end
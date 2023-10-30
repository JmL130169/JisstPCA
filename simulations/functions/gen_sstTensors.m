function [X, Y, u, V, W] = gen_sstTensors(p, q, N, rx, ry, K, orthogonal, Dx, Dy)
    u = cell(K, 1); V = cell(K, 1); W = cell(K, 1);
    X = tensor(zeros(p, p, N)); Y = tensor(zeros(q, q, N));
    if orthogonal
        V_mat = orth(normrnd(0, 1, [p, sum(rx)]));
        W_mat = orth(normrnd(0, 1, [q, sum(ry)]));
        u_mat = orth(normrnd(0, 1, [N 2]));
        start_V = 1; start_W = 1;
        for k = 1 : K
            u{k} = u_mat(:, k);
            V{k} = V_mat(:, start_V : (start_V + rx(k) - 1));
            W{k} = W_mat(:, start_W : (start_W + ry(k) - 1));
            start_V = start_V + rx(k); start_W = start_W + ry(k);
            X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(u{k})));
            Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(u{k})));
        end
    else
        for k = 1 : K
            u{k} = orth(normrnd(0, 1, [N 1]));
            V{k} = orth(normrnd(0, 1, [p, rx(k)]));
            W{k} = orth(normrnd(0, 1, [q, ry(k)]));
            X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(u{k})));
            Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(u{k})));
        end
    end
end
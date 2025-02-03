function [X, Y, ux, uy, V, W] = gen_sstTensors_heter(p, q, N, rx, ry, K, orthogonal, Dx, Dy, ratio)

    ux = cell(K, 1); ux2 = cell(K, 1); uy = cell(K, 1); V = cell(K, 1); W = cell(K, 1);
    X = tensor(zeros(p, p, N)); Y = tensor(zeros(q, q, N));

    if orthogonal
        V_mat = orth(normrnd(0, 1, [p, sum(rx)]));
        W_mat = orth(normrnd(0, 1, [q, sum(ry)]));
        ux_mat = orth(normrnd(0, 1, [N K]));
        uy_mat = orth(normrnd(0, 1, [N K]));
        start_V = 1; start_W = 1;
        for k = 1 : K
            ux{k} = ux_mat(:, k);
            uy{k} = ratio * ux_mat(:, k) + sqrt(1 - ratio^2) * uy_mat(:, k);
            uy{k} = uy{k} / norm(uy{k});
            V{k} = V_mat(:, start_V : (start_V + rx(k) - 1));
            W{k} = W_mat(:, start_W : (start_W + ry(k) - 1));
            start_V = start_V + rx(k); start_W = start_W + ry(k);
            X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(ux{k})));
            Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(uy{k})));
        end
    else
        for k = 1 : K
            ux{k} = orth(normrnd(0, 1, [N 1]));
            ux2{k} = orth(normrnd(0, 1, [N 1]));
            uy{k} = ratio * ux{k} + sqrt(1 - ratio^2) * ux2{k};
            uy{k} = uy{k} / norm(uy{k});
            V{k} = orth(normrnd(0, 1, [p, rx(k)]));
            W{k} = orth(normrnd(0, 1, [q, ry(k)]));
            X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(ux{k})));
            Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(uy{k})));
        end
    end

end
function [X, Y, u, V, W] = gen_sstTensors_K(p, q, N, rx, ry, Kx, Ky, orthogonal, Dx, Dy)

    V = cell(Kx, 1); W = cell(Ky, 1);
    X = tensor(zeros(p, p, N)); Y = tensor(zeros(q, q, N));
    K = max(Kx, Ky); u = cell(K, 1); 

    if orthogonal

        V_mat = orth(normrnd(0, 1, [p, sum(rx)]));
        W_mat = orth(normrnd(0, 1, [q, sum(ry)]));
        u_mat = orth(normrnd(0, 1, [N K]));
       
        start_V = 1; start_W = 1;
        for k = 1 : K
            u{k} = u_mat(:, k);
        end
        for k = 1 : Kx
            V{k} = V_mat(:, start_V : (start_V + rx(k) - 1));
            start_V = start_V + rx(k);
            X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(u{k})));
        end
        for k = 1 : Ky
            W{k} = W_mat(:, start_W : (start_W + ry(k) - 1));
            start_W = start_W + ry(k);
            Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(u{k})));
        end

    else

        for k = 1 : K
            u{k} = orth(normrnd(0, 1, [N 1]));
        end
        for k = 1 : Kx
            V{k} = orth(normrnd(0, 1, [p, rx(k)]));
            X = X + squeeze(ttt(tensor(V{k} * Dx{k} * V{k}'), tensor(u{k})));
        end
        for k = 1 : Ky
            W{k} = orth(normrnd(0, 1, [q, ry(k)]));
            Y = Y + squeeze(ttt(tensor(W{k} * Dy{k} * W{k}'), tensor(u{k})));
        end

    end

end
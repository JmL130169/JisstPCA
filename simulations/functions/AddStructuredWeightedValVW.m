function results = AddStructuredWeightedValVW(results, V, W, Dx, Dy, methodName)
    %% add results of estimated V, W to table "results"
    K = size(V, 1);
    p = size(V{1}, 1);
    q = size(W{1}, 1);
    for k = 1 : K
        method = repmat(methodName, p^2 + q^2, 1);
        FactorName = [repmat("V", p^2, 1); repmat("W", q^2, 1)];
        FactorNumber = repmat(k, p^2 + q^2, 1);
        RowInd = [reshape(repmat(1 : p, p, 1)', p^2, 1); reshape(repmat(1 : q, q, 1)', q^2, 1)];
        ColInd = [reshape(repmat(1 : p, p, 1), p^2, 1); reshape(repmat(1 : q, q, 1), q^2, 1)];
        Val = [reshape(V{k} * Dx{k} * V{k}', p^2, 1); reshape(W{k} * Dy{k} * W{k}', q^2, 1)];
        results = [results; table(method, FactorName, FactorNumber, RowInd, ColInd, Val)];
    end
end
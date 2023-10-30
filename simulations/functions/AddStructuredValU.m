function results = AddStructuredValU(results, u, methodName)
    %% add results of estimated V, W to table "results"
    K = size(u, 1);
    N = size(u{1}, 1);
    for k = 1 : K
        method = repmat(methodName, N, 1);
        FactorNumber = repmat(k, N, 1);
        RowInd = (1 : N)';
        Val = u{k};
        results = [results; table(method, FactorNumber, RowInd, Val)];
    end
end
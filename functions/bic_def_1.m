function [out1] = bic_def_1(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, deflation)
    
    [bic_rx, ~] = bic_def(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, deflation);
    out1 = bic_rx;
    
end

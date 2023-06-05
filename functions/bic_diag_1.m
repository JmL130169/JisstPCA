function [out1] = bic_diag_1(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, deflation)
    
    [bic_rx, ~] = bic_diag(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, deflation);
    out1 = bic_rx;
    
end

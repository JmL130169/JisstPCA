function [out2] = bic_def_2(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, deflation)
    
    [~, bic_ry] = bic_def(X, Y, rank_max, K, u0, lambda, tol, max_iter, method, deflation);
    out2 = bic_ry;
    
end

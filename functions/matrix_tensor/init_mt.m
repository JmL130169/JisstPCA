% This function outputs spectral initialization u0 when we input X and Y
% are tensor-matrix factorization

function [u0] = init_mt(X, Y)
    
    % matricization X by the third-mode and concatenate it with Y
    M = [double(tenmat(X, 3, [1, 2])), double(Y')];

    % using the leading singular vector as spectral initialization
    [a, ~, ~] = svds(M, 1);
    u0 = a;

end
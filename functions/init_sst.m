% This function outputs spectral initialization u0 when we input X

% Note: svd function in MATLAB automatically arrange the singular
% vectors/values in an decreasing order, we arrange that by hand to make
% sure we select the singular vectors with regard to larges singular values

function [u0] = init_sst(X)
    
    % matricization X by the third-mode and concatenate them
    M = double(tenmat(X, 3, [1, 2]));

    % using the leading singular vector as spectral initialization
    [a, ~, ~] = svds(M, 1);
    u0 = a;

end
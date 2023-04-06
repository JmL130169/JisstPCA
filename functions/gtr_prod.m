% This function is for trace product with regard to a tensor X, low-rank
% orthogonal matrix V and diagonal matrix D

% gtr_prod function is used in dJisst_single when signal is diagonal matrix
% D instead of scaler d

function u = gtr_prod(X, V, D)

    % The input of gtr_prod is p-p-N tensor X, p-r orthogonal matrix V and
    % p-p diagonal matrix D, the output of gtr_prod is a vector u of length N

    % For ith element of u, it is the trace of V'*X(:, :, i)*V*D, where
    % X(:, :, i) is the ith slice of tensor X
    sz = size(X, 3);
    u = zeros(sz, 1);
    for i = 1:sz
        u(i) = trace(V'*double(X(:, :, i))*V*D);
    end
end
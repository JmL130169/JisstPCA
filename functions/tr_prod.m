% This function is trace product used in Jisst_single

% input of tr_prod:
% X: semi-symmetric tensor of dimension p-p-N
% V: matrix of dimension p-r

% output of tr_prod:
% u: vector of length N

function u = tr_prod(X, V)
    % The input of tr_prod are a p-by-p-by-N tensor X and a p-by-k matrix
    % V, and the output of tr_prod is a N-dimensional vector, with the ith
    % element being the inner product of VV' and the ith slice of X (which
    % is X(:, :, i)).
    %
    % Example:
    % p = 20;
    % k = 5;
    % N = 50;
    % V = normrnd(0, 1, [p, k]);
    % V = orth(V);
    % u = normrnd(0, 1, [N, 1]);
    % u = u/(norm(u));
    % X = ttt(tensor(V*V'), tensor(u));
    % X = squeeze(X);
    % c = tr_prod(X, V); %<-- the output is a vector in R^{N} space.
    
    sz = size(X, 3);
    u = zeros(sz, 1);
    for i=1:sz
        u(i) = ttt(tensor(X(:,:,i)), tensor((V*V')), [1:2]);
    end
end

% sin theta distance using operator norm

% input of sin_do:
% X: low-rank matrix of dimension p-rx
% Y: low-rank matrix of dimension p-ry

% output of sin_do:
% d: sin_theta distance between the subspace spanned by X and Y, which is operator norm of sin_theta diagonal matrix

function d = sin_do(X, Y)
    [~, D, ~] = svd(X'*Y);
    D = diag(D);
    ac = acos(D);
    s = sin(ac);
    d = norm(diag(s));
end

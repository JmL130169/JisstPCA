% sin theta distance using operator norm

function d = sin_do(X, Y)
    [~, D, ~] = svd(X'*Y);
    D = diag(D);
    ac = acos(D);
    s = sin(ac);
    d = norm(diag(s));
end
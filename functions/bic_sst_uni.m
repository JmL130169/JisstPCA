% This function calculate value of BIC when the input rank is rx in
% single factor case. We always select rx that have smallest BIC as the
% input ranks.

% input of bic_sst_uni:
% rx: possible rank of single factor X
% X: single factor semi-symmetric tensor
% u0: initialization for implementing sst_single
% tol, max_iter: inputs that needed in sst_single

% output of bic_sst_uni:
% bic_val: the corresponding BIC value for rank rx, given data tensor X

function [bic_val] = bic_sst_uni(rx, X, u0, tol, max_iter)
    
    sz_x = size(X);

    n = sz_x(1)*sz_x(2)*sz_x(3);
   
    pe = sum(sz_x(1)*rx + sz_x(3)); % number of efficient parameters
    
    [~, ~, ~, hat_X] = sst_single(X, u0, rx, tol, max_iter);

    % log likelihood function of X based on estimation using rank rx
    Tx = n * log(mean(double((X - hat_X).^2), 'all'));

    bic_val = pe*log(n) + Tx;

end
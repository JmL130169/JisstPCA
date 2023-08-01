% This function calculate value of BIC when the input rank is rx and ry in
% single factor case. We always select rx, ry that have smallest BIC as the
% input ranks.

% input of bic_sst:
% rx, ry: possible rank of single factor X and Y
% X, Y: single factor semi-symmetric tensors
% u0: initialization for implementing Jisst_single
% lambda: scaler, prespecified or data-driven
% tol, max_iter: inputs that needed in Jisst_single

% output of bic_sst:
% bic_val: the corresponding BIC value for rank rx and ry, given data X and Y

function [bic_val] = bic_dsst(rx, ry, X, Y, u0, lambda, tol, max_iter)
    
    sz_x = size(X);
    sz_y = size(Y);

    n1 = sz_x(1)*sz_x(2)*sz_x(3);
    n2 = sz_y(1)*sz_y(2)*sz_y(3);
    n = n1 + n2; % number of observation data

    pe = sum(sz_x(1)*rx + sz_y(1)*ry + sz_x(3)); % number of efficient parameters
    
    [~, ~, ~, ~, ~, hat_X, hat_Y] = dJisst_single(X, Y, u0, rx, ry, lambda, tol, max_iter);

    % log likelihood function of X based on estimation using rank rx
    Tx = n1 * log(mean(double((X - hat_X).^2), 'all'));

    % log likelihood function of y based on estimation using rank ry
    Ty = n2 * log(mean(double((Y - hat_Y).^2), 'all'));

    bic_val = pe*log(n) + Tx + Ty;
end
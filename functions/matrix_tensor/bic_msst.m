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

function [bic_val] = bic_msst(r, X, Y, u0, lambda, tol, max_iter)
    
    sz_x = size(X);

    n = sz_x(1)*sz_x(2)*sz_x(3); % number of observation data

    pe = sum(sz_x(1)*r + sz_x(3)); % number of efficient parameters
    
    [~, ~, ~, hat_X, ~, ~, ~] = Jimt_single(X, Y, u0, r, lambda, tol, max_iter);

    % log likelihood function of X based on estimation using rank rx
    Tx = (-0.5)*(log(2*pi) + (X - hat_X).^2);

    bic_val = pe*log(n) - 2*sum(double(Tx), 'all');
end
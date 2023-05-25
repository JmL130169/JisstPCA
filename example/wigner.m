function w = wigner(n,isreal)
%WIGNER The Wigner matrix.
% WIGNER(N,ISREAL) generates an N x N symettric Wigner matrix.
% If ISREAL = 1 then the elements of W are real.
% If ISREAL = 0 then the elements of W are complex. %
% N is an integer while ISREAL equals either 0 or 1.
% WisanNxNmatrix. %
% References:
% [1] Alan Edelman, %
%
% [2] Alan Edelman,
% [3] E. P. Wigner, %
%
%
% Alan Edelman and Raj Rao, Sept. 2004.
% $Revision: 1.0 $ $Date: 2004/09/10 23:21:18 $
if(isreal==1)
    g = randn(n,n);
else
    g = (randn(n,n) + i*randn(n,n))/sqrt(2);
end
    w = (g + g') / sqrt(2);
end
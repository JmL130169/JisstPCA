function w = wigner_sst(p, N)
%semi-symmetric tensor with each slice generated from symmetric Gaussian
%ensemble noise
g = randn(p, p, N);
w = (g + permute(g, [2 1 3])) / sqrt(2);
end
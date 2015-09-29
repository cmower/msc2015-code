function [G, T] = rand_acor_2b(n, k)
%RAND_ACOR_2B Random approximate correlation matrix, two block variant.
%   [G, T] = RAND_ACOR_2B(n, k) produces an n-by-n approximate correlation
%   matrix of the form G = [A Y; Y' B] where A and B are valid correlation
%   matrices. k is the size of the matrix A. The matrix T = blkdiag(A, B)
%   is the target in the shrinking method.
%   k defaults to ceil(n/2).
%
%   By C. E. Mower, 26/07/2015.
%

% setup
if nargin < 2, k = ceil(n/2); end

A = gallery('randcorr', k);
B = gallery('randcorr', n-k);

while 1
    Y = -1 + 2*randomNumbers(k, n-k, 'u');
    G = [A Y; Y' B];
    [~, info] = chol(G);
    if info, break; end
end

T = blkdiag(A, B);

end
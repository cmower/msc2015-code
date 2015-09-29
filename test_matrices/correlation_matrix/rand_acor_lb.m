function [G, T] = rand_acor_lb(n, k, r)
%RAND_ACOR_LB Random Approximate Correlation Matrix, leading block variant.
%   [G, T] = RAND_ACOR_LB(n, k, r) produces an n-by-n random approximate
%   correlation matrix with a k-by-k positive (semi)definite leading block.
%   The matrix G has the form [A Y;Y' B] where A is the k-by-k valid 
%   correlation matrix of rank r. The matrix T = blkdiag(A, I) is the
%   target matrix used in the shrinking method, I is the order n-k identity
%   matrix.
%   k defaults to ceil(n/2).
%   r defaults to k.
%
%   See also gallery.
%
%   By C. E. Mower, 26/07/2015.
%

% setup
if nargin < 2, k = ceil(n/2); end
if nargin < 3, r = k; end

% generate evals of A
evals_zero = zeros(1, k-r);
while 1
    evals_mid = rand(1, r - 1);
    last_eval = k - sum(evals_mid);
    if last_eval > 0, break; end
end
evals = [evals_zero evals_mid last_eval];

A = gallery('randcorr', evals);

% generate G
while 1
    Y = -1 + 2*randomNumbers(k, n-k, 'u');
    B = rand_acor(n-k);
    G = [A Y; Y' B];
    [~, info] = chol(G);
    if info, break; end
end

T = blkdiag(A, eye(n-k));

end
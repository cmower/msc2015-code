function G = rand_acor(n)
%RAND_ACOR Random approximate correlation matrix.
%   G = RAND_ACOR(n) produces an n-by-n symmetric indefinite approximate
%   correlation matrix. Note, n must be greater than 2.
%
%   See also rand_sym.
%
%   By C. E. Mower, 26/07/2015.
%

while 1
    G = -1 + 2*rand_sym(n, 'u');
    G(1:n+1:n^2) = 1; % unit diagonal
    [~, info] = chol(G);
    if info, break; end
end

end

function G = rand_sym(x, p)
%RAND_SYM Random Symmetric Matrix.
%   G = RAND_SYM(n, p) produces a random symmetric matrix G with elements
%   from a distribution specified by p. 
%   p = 'n' elements of X are from a normal distribution (default).
%   p = 'u' elements of X are from a uniform distribution, i.e. from (0,1).
%   G = RAND_SYM(x) produces a random symmetric matrix having eigenvalues
%   given by the vector x, where length(x) > 1.
%
%   See also randomNumbers, randn, rand.
%
%   By C. E. Mower, 03/08/2015.
%

if nargin < 2, p = 'n'; end

n = length(x);

if n == 1
    G = randomNumbers(x, x, p);
else
    D = diag(x);
    X = rand_sym(n, p);
    Q = orth(X);
    G = Q*D*Q';
end

G = (G + G') / 2;

end
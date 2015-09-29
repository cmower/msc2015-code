function G = TD_Toep(n, gamma)
%TD_TOEP Symmetric tridiagonal Toeplitz Matrix with unit diagonal.
%   G = TD_TOEP(n) is the matrix with diagonal 1, and both super and sub 
%   diagonal 1.
%   G = TD_TOEP(n, gamma) is the tridiagonal Toeplitz matrix G of 
%   order n with diagonal 1, and both super and sub diagonals gamma.
%
%   By C. E. Mower, 17/7/2015.
%

if nargin < 2, gamma = 1; end

G = full(gallery('tridiag', n, gamma, 1, gamma));

end
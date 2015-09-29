function r = iscorr(A)
%ISCORR True if A is a correlation matrix.
%   ISCORR(A) returns true (logical 1) if A is a correlation matrix and 
%   false (logical 0) otherwise.
%
%   NOTE:
%   * Correlation matrices are defined to be positive semidefinite (psd), 
%     this function tests definiteness by a Cholesky Factorization. If you
%     believe your matrix to be psd then this will likely fail. 
%     Solution: compute the minimum eigenvalue using min(eig(A)) and check
%     if it is equal to zero (or within some tolerance).
%
%   By C. E. Mower, 19/10/2014.
%

n = length(A);
tol = n*eps;

% check range of elements.
maxA = max(max(A));
minA = min(min(A));
r = minA > -1 && maxA < 1;
if ~r, 
    r = abs(minA + 1) > -tol;
    if ~r, return; end
    r = abs(maxA - 1) > -tol;
    if ~r, return; end
end

% check unit diagonal
for i = 1:n
    d = abs(A(i,i) - 1);
    r = d < n*eps;
    if ~r, return; end
end

% check for pd (and symmetry).
r = ispd(A);

end


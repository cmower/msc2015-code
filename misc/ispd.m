function r = ispd(A)
%ISPD True if A is symmetric positive definite.
%   ISPD(A) returns true (logical 1) if A is symmetric positive definite 
%   and false (logical 0) otherwise.
%
%   By C. E. Mower, 09/08/2015.
%

if ~issymmetric(A), error('A must be symmetric.'); end
[~, p] = chol(A);
r = p == 0;

end


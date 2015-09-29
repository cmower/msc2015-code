function r = ispsd(A)
%ISPSD True if A is symmetric positive semidefinite.
%   ISPSD(A) returns true (logical 1) if A is symmetric positive 
%   semidefinite and false (logical 0) otherwise.
%
%   See also, eps.
%
%   By C. E. Mower, 09/08/2015.
%

if ~issymmetric(A), error('A must be symmetric.'); end
r = min(eig(A)) > -eps;
    
end


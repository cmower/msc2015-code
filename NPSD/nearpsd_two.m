function X = nearpsd_two(G)
%NEARPSD_TWO Nearest Positive Semidefinite Matrix in the 2 norm.
%   X = NEARPSD_TWO(G) computes the nearest positive semidefinite matrix in
%   the 2 norm to the symmetric indefinite matrix G. See [1].
%
%   References:
%   [1] Nicholas J. Higham. Computing a nearest symmetric positive 
%       semidefinite matrix. Linear Algebra Appl., 103:103?118, 1988.
%
%   By C. E. Mower, 03/08/2015.
%

emin = min(eig(G)); % note for indefinite G this is negative.
X = G - emin*eye(length(G));

end


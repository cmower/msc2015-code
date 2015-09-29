function X = nearpsd_fro(G)
%NEARPSD_FRO Nearest Positive Semidefinite Matrix in the Frobenius norm.
%   X = NEARPSD_FRO(G) computes the nearest positive semidefinite matrix in
%   the frobenius norm to the symmetric indefinite matrix G. See [1].
%
%   References:
%   [1] Nicholas J. Higham. Computing a nearest symmetric positive 
%       semidefinite matrix. Linear Algebra Appl., 103:103-118, 1988.
%
%   By C. E. Mower, 03/08/2015.
%

[Q, D] = eig(G);
X = Q*diag(max(diag(D),0))*Q';
X = (X + X') / 2; % ensure symmetry

end


function A = cor_bar(P)
%COR_BAR  Calculates approximate sample correlation matrix.
%   A=COR_BAR(P) Produces an n-by-n approx correlation matrix based on
%   data of size m-by-n. 
%   n columns of different random variables observed at m different times.
%   P has missing data represented by NaNs.
%
%   References:
%       [1] Craig Lucas, Computing the nearest covariance and correlation
%           matrices, M.Sc. Thesis, University of Manchester, Manchester,
%           England, October 2001. 68pp.
%   
%   By C. Lucas 2001.
%
%   Modified by C. E. Mower, 01/04/2015.
%
%   See also, cov_bar.
%

[~,n] = size(P);
S=cov_bar(P);
D=diag(1./sqrt(diag(S)));
A=D*S*D;
A = (A + A') / 2; A(1:n+1:n^2) = 1; % ensure symmetry and unit diagonal.
end
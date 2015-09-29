function S = cov_bar(P)
%COV_BAR  Calculates approximate sample covariance matrix.
%   S = COV_BAR(P) Produces an n-by-n approx covariance matrix based on
%   data of size m-by-n. 
%   n columns of different random variables observed at m different times.
%   P has missing data represented by NaNs.
%
%   References:
%       [1] Craig Lucas, Computing the nearest covariance and correlation
%           matrices, M.Sc. Thesis, University of Manchester, Manchester,
%           England, October 2001. 68pp.
%
%   By C. Lucas, 2001.
%
%   Modified by C. E. Mower, 01/04/2015.
%

[~,n] = size(P);

S = zeros(n);

for i = 1:n
    xi = P(:,i);

    for j=1:i
    xj = P(:,j);

        % create mask for data values that are 'common'
        p = ~isnan(xi) & ~isnan(xj);

        S(i,j) = (xi(p) - mean(xi(p)))'*(xj(p) - mean(xj(p)));

        % normalise over effective sample size i.e. sum(p)-1
        S(i,j) = 1/(sum(p)-1)*S(i,j);

        S(j,i) = S(i,j);

    end
end


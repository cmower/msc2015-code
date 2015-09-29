function [alpha, iter] = shrink_bisect_Hme(G, H, theta, tol)
%SHRINK_BISECT_HME Shrinking via bisection; elementwise weighting and
%   optionally specifies a minimum eigenvalue.
%   [alpha, iter] = SHRINK_BISECT_HME(G, H, theta, tol) computes the
%   smallest alpha in [0, 1] such that S = alpha*(H.*G) + (1 - alpha)G is 
%   positive (semi)definite using the bisection variant of the shrinking 
%   method. Elements may be fixed, or weighted, by specifying them in the 
%   matrix H: if it is required that S(i,j) = G(i,j) then specify
%   H(i,j) = 1. If the element S(i,j) is to be weighted then specify 
%   H(i,j) = h, where h is a number from the open interval (0,1).
%   Otherwise, if no weighting is required, set H(i,j) = 0.
%   Choosing theta from the open interval (0,1) enforces positive
%   definiteneess of S by setting a lower bound on the minimum eigenvalue,
%   denoted phi. Note that phi = theta*min(eig(T)), where T = H.*G.
%   Choosing theta = 0 computes an S that is positive semidefinite.
%
%   NOTE:
%   * There are NO DEFAULTS.
%   * Input is NOT check for validity.
%
%   By C. E. Mower, 30/07/2015.
%

% setup for bisection
T = H.*G;
I = eye(size(G));
phi = theta*min(eig(T));
left = 0;
right = 1;
err = 1;
iter = 0;

while err > tol
    
    iter = iter + 1;
    
    mid = (left + right) / 2;
    
    S = mid*T + (1 - mid)*G - phi*I;
    
    [~, p] = chol(S);
    
    if p
        left = mid;  % S is indefinite.
    else
        right = mid; % S is pos-def
    end
    
    err = right - left;
    
end

alpha = right;

end


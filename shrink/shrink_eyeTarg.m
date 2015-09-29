function alpha = shrink_eyeTarg(G, theta)
%SHRINK_EYETARG Shrinking with identity target.
%   alpha = SHRINK_EYETARG(G, theta) computes the smallest alpha in [0,1]
%   such that S(alpha) = alpha*I + (1 - alpha)*G is positive
%   (semi)definite, where G is symmetric indefinite and T is symmetric
%   positive (semi)definite. Choosing theta from the open interval (0,1) 
%   enforces positive definiteneess of S by setting a lower bound on the 
%   minimum eigenvalue. Choosing theta = 0 computes an S that is positive 
%   semidefinite.
%
%   NOTE:
%   * There are NO DEFAULTS.
%   * Input is NOT checked for validity.
%   
%   By C. E. Mower, 06/07/2015.
%

emin = min(eig(G));
alpha = (theta - emin) / (1 - emin);

end


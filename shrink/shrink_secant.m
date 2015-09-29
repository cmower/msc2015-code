function [alpha, iter] = shrink_secant(G, T, maxit, tol)
%SHRINK_SECANT Shrinking via Sectant method.
%   [alpha, iter] = SHRINK_SECTANT(A, T, tol) computes the smallest alpha
%   in [0, 1] such that S(alpha) = alpha*T + (1 - alpha)*G is positive
%   semidefinite, where G is symmetric indefinite and T is symmetric
%   positive (semi)definite. iter is the number of iterations of the
%   Secant method, maxit is the maximum number of iterations performed and
%   tol is a convergence tolerance.
%
%   NOTE:
%   * There are NO DEFAULTS.
%   * Input is NOT check for validity.
%
%   By C. E. Mower, 30/06/2015.
%

% intial values.

alpha = zeros(1, maxit);
f = zeros(1, maxit);
alpha(1) = 0; f(1) = min(eig(G));
alpha(2) = 1; f(2) = min(eig(T));


for iter = 3:maxit
    
    alpha(iter) = alpha(iter-1) - ...
        (f(iter - 1))*((alpha(iter-1) - alpha(iter-2)) / ...
                                                  (f(iter-1) - f(iter-2)));
    
    err = abs(alpha(iter) - alpha(iter - 1));
    
    if err < tol, break; end
    
    S = alpha(iter)*T + (1 - alpha(iter))*G;
    f(iter) = min(eig(S));
end

alpha = alpha(iter);

end


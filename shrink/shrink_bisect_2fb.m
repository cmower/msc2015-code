function [alpha, iter] = shrink_bisect_2fb(G, k, tol)
%SHRINK_BISECT_2FB Shrinking by bisection method; two fixed block variant.
%   [alpha, iter] = SHRINK_BISECT_2FB(G, k, tol) computes the smallest 
%   alpha in [0,1] such that alpha*T + (1 - alpha)*G is positive 
%   semidefinite, where G = [A Y; Y' B] is the order n symmetric 
%   indefinite matrix, A is a k-by-k positive definite leading block, B is 
%   an (n-k)-by-(n-k) positive definite matrix and Y is a k-by-(n-k) full 
%   matrix. iter is the number of iterations of bisection and tol is a 
%   convergence tolerance.
%
%   NOTE:
%   * There are NO DEFAULTS.
%   * Input is NOT check for validity.
%
%   By C. E. Mower, 30/6/2015.
%

% setup for for bisection
left = 0;
right = 1;
err = right - left;
iter = 0;

% setup for bisection
n = length(G);
R11 = chol(G(1:k,1:k));
X = R11'\G(1:k,k+1:n);
X = X'*X;

% bisection
while err > tol
    
    iter = iter + 1;
    
    % build P
    mid = 0.5*(left + right);
    P = G(k+1:n, k+1:n) - (1-mid)^2*X;
    
    % attempt cholesky factorization
    [~,p] = chol(P);

    if p 
       left = mid; % S is indefinite
    else
       right = mid;
    end

    err = right - left;
    
end

alpha = right;

end

    


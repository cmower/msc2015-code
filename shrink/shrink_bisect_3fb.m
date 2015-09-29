function [alpha, iter] = shrink_bisect_3fb(G, k1, k2, tol)
%SHRINK_BISECT_3FB Shrinking via bisection method; 3 fixed block variant.
%   [alpha, iter] = SHRINK_BISECT_3FB(A, k1, k2, tol) computes the smallest 
%   alpha in [0,1] such that S(alpha) = alpha*T + (1 - alpha)*G is positive 
%   semidefinite, where G is symmetric indefinite with unit diagonal and 
%   T = blkdiag(A1, A2, A3) is the target matrix. See [1, Sec. 2.1.9].
%   iter is the number of iterations of bisection and tol is a convergence 
%   tolerance.
%
%   NOTE:
%   * There are NO DEFAULTS.
%   * Input is NOT check for validity.
%
%   By C. E. Mower, 01/04/2015.
%

% Slice and dice
n = length(G);
A1 = G(1:k1, 1:k1);
A2 = G(k1+1:k1+k2, k1+1:k1+k2);
A3 = G(k1+k2+1:n, k1+k2+1:n);

Y1 = G(1:k1,k1+1:k1+k2);
Y2 = G(1:k1,k1+k2+1:n);
Y3 = G(k1+1:k1+k2,k1+k2+1:n);

% setup for bisection
left = 0;
right = 1;
err = 1;
iter = 0;
 
R11 = chol(A1);
X1 = R11' \ Y1;
X2 = R11' \ Y2;

Z1 = X1'*X1;
Z2 = X2'*X2;
Z3 = X1'*X2;

% bisection
while err > tol
    
    iter = iter + 1;
    
    mid = 0.5*(left + right);

    P1 = A2 - (1 - mid)^2*Z1;
    
    [R22, ifail] = chol(P1);
    
    if ifail
        left = mid;
    else
        R23 = R22' \ ( (1 - mid)*Y3 - (1 - mid)^2*Z3);
        P2 = A3 - (1 - mid)^2*Z2 - R23'*R23;
        [~, ifail] = chol(P2);
        if ifail
            left = mid;
        else
            right = mid;
        end
    end

    err = right - left;
    
end

alpha = right;

end

    


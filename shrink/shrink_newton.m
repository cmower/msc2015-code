function [alpha, iter] = shrink_newton(G, T, maxit, tol)
%SHRINK_NEWTON Shrinking via Newton's method.
%   [alpha, iter] = SHRINK_NEWTON(G, T, maxit, tol) computes the smallest
%   alpha in [0, 1] such that S(alpha) = alpha*T + (1 - alpha)*G is
%   positive semidefinite, where G is symmetric indefinite and T is
%   symmetric positive (semi)definite. iter is the number of iterations of
%   the Newton method and tol is a convergence tolerance.
%
%   NOTE:
%   * There are NO DEFAULTS.
%   * Input is NOT checked for validity.
%
%   By N. J. Higham and N. Strabic, 2014.
%
%   Modified by C. E. Mower, 06/04/2015.
%

% Copyright (c) 2014, Nicholas J. Higham and Natasa Strabic.
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% * Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 
% * Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

% setup for Newton
x0 = 0; % Starting point
M = G - T;
iter = 0;

% Newton
for iter = 0:maxit
    S = x0*T + (1-x0)*G;
    [V, D] = eig(S);
    x_min = V(:,1); % Eigenvector for the smallest eigenvalue
    x1 = x_min'*G*x_min/(x_min'*M*x_min);
    if abs(x0-x1) <= tol, alpha = x1; return, end
    x0 = x1;
    if iter == maxit
        error('Not converged in %2.0f iterations', maxit)
    end
end

end
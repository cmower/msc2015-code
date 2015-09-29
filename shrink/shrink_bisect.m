function [alpha, iter] = shrink_bisect(G, T, tol)
%SHRINK_BISECT Shrinking via the bisection method.
%   [alpha, iter] = SHRINK_BISECT(G, T, tol) computes the smallest alpha
%   in [0,1] such that alpha*T + (1 - alpha)*G is positive semidefinite, 
%   where G is symmetric indefinite and T is a positive definite target 
%   matrix. iter is the number of iterations of bisection and tol is a 
%   convergence tolerance.
%
%   NOTE:
%   * There are NO DEFAULTS.
%   * Input is NOT check for validity.
%
%   By N. J. Higham and N. Strabic, 2014.
%
%   Modified by C. E. Mower, 01/04/2015.
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

% setup for bisection
left = 0;
right = 1; 
err = 1;
iter = 0;

% bisection
while err > tol
    
   iter = iter + 1;

   middle = 0.5*(left + right);
   
   % build S  
   S = middle*T + (1 - middle)*G;

   % attempt cholesky factorization
   [~, p] = chol(S);
   
   if p
      left = middle;  % S is indefinite.
   else
      right = middle; % S is pos-def
   end

   err = right - left;

end

alpha = right;

end
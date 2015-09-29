function [alpha, iter] = shrink_bisect_fb(G, k, tol)
%SHRINK_BISECT_FB Shrinking by bisection method; fixed block variant.
%   [alpha, iter] = SHRINK_BISECT_FB(G, k, tol) computes the smallest alpha
%   in [0,1] such that alpha*T + (1 - alpha)*G is positive semidefinite, 
%   where G = [A Y; Y' B] is symmetric indefinite, A is k-by-k positive
%   definite and T is the target matrix. tol is a convergence tolerance.
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

% Slice and dice
A = G(1:k, 1:k);
B = G(k+1:end, k+1:end);
Y = G(1:k, k+1:end);
I = eye(size(B));

left = 0;
right = 1;
err = right - left;
iter = 0;
 
R11 = chol(A);
X = R11'\Y;
Z = X'*X;

while err > tol
    
    iter = iter + 1;
    
    mid = 0.5*(left + right);
    P = mid*I + (1-mid)*B - (1-mid)^2*Z;

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

    


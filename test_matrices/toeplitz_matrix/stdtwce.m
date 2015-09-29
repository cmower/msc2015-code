function G = stdtwce(n, gamma, beta)
%STDTWCE Symmetric Tridiagonal Toeplitz Matrix with unit diagonal and 
%   corner elements.
%   G = stdtwce(n, gamma, beta) is the order n symmetric matrix with ones 
%   along the main diagonal, gamma super and sub diagonal and corner
%   elements beta.
%
%   By C. E. Mower, 10/08/2015.

if n == 1, G = 1; return; end
if n == 2, G = [1 gamma; gamma 1]; return; end

G = TD_Toep(n, gamma);
G(1,n) = beta;
G(n,1) = beta;

end


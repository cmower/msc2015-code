function G = KMS_Toep(n, rho)
%KMS_TOEP Kac-Murdock-Szego Toeplitz Matrix.
%   G = KMS_TOEP(n, rho) returns the KMS Toeplitz matrix G of order n where
%   G(i,j) = rho^abs(i - j) and rho is a real number. 
%   
%   By C. E. Mower, 10/08/2015.
%

G = gallery('kms', n, rho);

end
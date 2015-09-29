function X = randomNumbers(m, n, p)
%RANDOMNUMBERS Matrix of Random Numbers.
%   X = RANDOMNUMBERS(m, n, p) produces a random m-by-n matrix with 
%   elements from a distribution specified by p. 
%   p = 'n' elements of X are from a normal distribution.
%   p = 'u' elements of X are from a uniform distribution, i.e. from (0,1).
%
%   See also randn, rand.
%
%   By C. E. Mower, 03/08/2015.
%

switch p
    case 'n'
        X = randn(m, n);
    case 'u'
        X = rand(m, n);
end

end


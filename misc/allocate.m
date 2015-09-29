function P = allocate(P, k, t)
%ALLOCATE Randomly allocates elements to a matrix.
%   P = ALLOCATE(P, k, t) allocates N t's to the matrix P. If 0 < k < 1
%   then N = ceil(m*n*k) and if k is a positive integer then N = k.
%
%   NOTE:
%   * k defaults if unspecified or k = [].
%   * t defaults to NaN.
%
%   By C. E. Mower, 02/04/2015.
%

% setup
if nargin < 2 || isempty(k), k = 0.3; end
if nargin < 3, t = NaN; end

% error check
if k < 0, error('k must either be a positive integer or from the set (0,1).'); end

% generate P
[m, n] = size(P);
if k < 1
    num_t = ceil(m*n*k);
else
    num_t = k;
end

while num_t > 0
    i = randi(m);
    j = randi(n);
    if P(i,j) ~= t
        P(i,j) = t;
        num_t = num_t - 1;
    end
end

end
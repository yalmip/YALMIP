function k = findrows(A, b)
%FINDROWS Find indices of a given row within a matrix.
%
%   FINDROWS(A, B) returns a column vector with the indices of the rows
%   in the matrix A that are identical to the row vector B.  If no rows
%   in A are identical to B, an empty vector is returned.
%
%   The methods uses a for-loop, but it uses less memory and is in many
%   cases a lot faster than the vectorized methods
%
%      find( all( A == repmat(b, size(A, 1), 1), 2 ) )
%      find( all( A == b(ones(size(A, 1), 1),:), 2 ) )
%
%   See also FIND, FINDCOLS.

%   Author:      Peter J. Acklam
%   Time-stamp:  2002-03-03 13:51:15 +0100
%   E-mail:      pjacklam@online.no
%   URL:         http://home.online.no/~pjacklam

if isempty(A)
    k = [];
    return
end

k = find( A(:,1) == b(1));
top = size(A, 2);
for j = 2:top
    k = k(A(k,j) == b(j));
    if isempty(k)
        return
    end
end
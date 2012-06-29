function index = end(X,position,numindices)
%END (overloaded)

% Author Johan Löfberg
% $Id: end.m,v 1.1 2006-07-13 19:40:59 joloef Exp $

if nargin == 2
    index = X.dim(position);
else
    if numindices == 1
        index = prod(X.dim);
    elseif numindices == length(X.dim)
        index = X.dim(position);
    else
        error('Indexing logic not supported yet');
    end
end
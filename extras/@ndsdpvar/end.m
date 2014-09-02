function index = end(X,position,numindices)
%END (overloaded)

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
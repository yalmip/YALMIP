function index = end(X,position,numindices)
%END (overloaded)

switch numindices
    case 1
        % User has written someting like X(end)
        sizes = X.dim;%[X.n X.m];
        if min(sizes)>1
            index = prod(sizes);
        else
            index = max(sizes);
        end
    case 2
        sizes = X.dim;%[X.n X.m];
        index = sizes(position);
    otherwise
        error('Indexation dimension cannot exceed 2');
end

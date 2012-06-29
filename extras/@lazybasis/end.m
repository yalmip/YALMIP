function index = end(X,position,numindices)
%END (overloaded)

% Author Johan Löfberg
% $Id: end.m,v 1.1 2005-10-12 16:05:54 joloef Exp $

switch numindices
    case 1
        % User has written someting like X(end)
        sizes = size(double(X));
        if min(sizes)>1
            index = prod(sizes);
        else
            index = max(sizes);
        end
    case 2
        sizes = size(double(X));
        index = sizes(position);
    otherwise
        error('Indexation dimension cannot exceed 2');
end

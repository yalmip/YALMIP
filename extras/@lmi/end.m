function index = end(X,position,numindices)
%end               Overloaded

switch numindices
    case 1
        % User has written someting like X(end)
        index = length(X);
    case 2
        index = length(X);
    otherwise
        error('Indexation dimension cannot exceed 2');
end
  
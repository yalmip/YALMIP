function index = end(X,position,numindices)
%end               Overloaded

% Author Johan Löfberg 
% $Id: end.m,v 1.1 2004-06-17 08:40:04 johanl Exp $   

  switch numindices
   case 1
    % User has written someting like X(end)    
    index = length(X);
   case 2
    index = length(X);
   otherwise
    error('Indexation dimension cannot exceed 2');
  end
  
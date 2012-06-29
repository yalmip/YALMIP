function F = ismember(x,F)
% Internal class for constraint list

% Author Johan Löfberg
% $Id: ismember.m,v 1.1 2004-06-17 08:40:03 johanl Exp $

F = replace(F,recover(depends(F)),x);

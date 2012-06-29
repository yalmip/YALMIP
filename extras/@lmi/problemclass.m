function p = problemclass(F,h)
% PROBLEMCLASS Derives an optimization object and determines the class
%

% Author Johan Löfberg
% $Id: problemclass.m,v 1.2 2004-07-19 13:54:35 johanl Exp $

if nargin < 2
    h = [];
end
[aux1,aux2,aux3,model] = export(F,h);
p = problemclass(model);

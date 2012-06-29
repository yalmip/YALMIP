function F = tag(F,text)
% TAG
%
% Sets the tag on a constraint.

% Author Johan Löfberg
% $Id: tag.m,v 1.1 2006-05-16 12:32:36 joloef Exp $

if nargin == 1
    F = F.clauses{1}.handle;
else
for i = 1:length(F.clauses)
    F.clauses{i}.handle = text;
end
end

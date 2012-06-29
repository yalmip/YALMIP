function F = colon(F,tag)
% COLON Overloaded

% Allows the syntax (x>0):Tag in order to give names/descriptions to
% constraints

% Author Johan Löfberg
% $Id: colon.m,v 1.1 2009-05-04 09:15:22 joloef Exp $

for i = 1:length(F.clauses)
    F.clauses{i}.handle = tag;
end
	
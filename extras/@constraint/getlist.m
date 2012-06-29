function [F,strict,LMIIdentifiers,tags] = getlist(X)
% Internal class for constraint lists

% Author Johan Löfberg
% $Id: getlist.m,v 1.5 2009-04-29 07:48:12 joloef Exp $

superiorto('sdpvar');
superiorto('double');

F = X.Evaluated;
strict = X.strict;
LMIIdentifiers = X.ConstraintID;
tags = X.tag;

% FIX : treat equalities better
if isequal(X.List{2},'==') 
    F{1}=sethackflag(F{1},3);
end

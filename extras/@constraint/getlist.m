function [F,strict,LMIIdentifiers,tags] = getlist(X)
% Internal class for constraint lists

F = X.Evaluated;
strict = X.strict;
LMIIdentifiers = X.ConstraintID;
tags = X.tag;

% FIX : treat equalities better
if isequal(X.List{2},'==') 
    F{1}=sethackflag(F{1},3);
end
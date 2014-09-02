function F = ismember(x,F)
% Internal class for constraint list

F = lmi(F);
F = ismember(x,F);


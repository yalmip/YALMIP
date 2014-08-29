function F = replace(F,x,w)
% Internal class for constraint list

F = lmi(F);
F = replace(F,x,w);


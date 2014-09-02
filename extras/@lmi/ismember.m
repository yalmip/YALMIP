function F = ismember(x,F)
% Internal class for constraint list

F = replace(F,recover(depends(F)),x);

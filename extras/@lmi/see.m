function see(F)
%see               Displays internal structure of matrix variable in a constraint

F = flatten(F);
if length(F.LMIid)>1
    disp('Can only apply SEE on set objects with one constraint');
else    
    see(F.clauses{1}.data);
end
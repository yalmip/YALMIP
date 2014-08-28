function itslinear=islinear(F)
%ISLINEAR Check if all constraints are linear       
%
% p = islinear(F)
%
% F : SET object
% p : boolean 0/1

[monomtable,variabletype] = yalmip('monomtable');
itslinear = 1;
i = 1;
F = flatten(F);
while itslinear && (i<=length(F.LMIid))
    Fi = F.clauses{i};
    xvars = getvariables(Fi.data);
    itslinear = itslinear & ~any(variabletype(xvars));
    i = i + 1;    
end
function itslinear=islinear(F)
%ISLINEAR Check if all constraints are linear       
%
% p = islinear(F)
%
% F : SET object
% p : boolean 0/1

% Author Johan Löfberg 
% $Id: islinear.m,v 1.3 2005-02-04 10:10:26 johanl Exp $  

[monomtable,variabletype] = yalmip('monomtable');
itslinear = 1;
i = 1;
while itslinear & (i<=length(F.clauses))
    Fi = F.clauses{i};
    xvars = getvariables(Fi.data);
    itslinear = itslinear & ~any(variabletype(xvars));
    i = i + 1;    
end
function itssigmonial=issigmonial(F)
%ISLINEAR Check if all constraints are linear       
%
% p = islinear(F)
%
% F : SET object
% p : boolean 0/1

% Author Johan Löfberg 
% $Id: issigmonial.m,v 1.2 2005-02-04 10:10:27 johanl Exp $  


monomtable = yalmip('monomtable');
xvars = [];

itssigmonial = 1;
i = 1;
while itssigmonial & (i<=length(F.clauses))
    Fi = F.clauses{i};
    xvars = getvariables(Fi.data);
    monomtableX = monomtable(xvars,:);
    YESNO = any(find(any(0>monomtableX,2) | any(monomtableX-fix(monomtableX),2)));   
    itssigmonial = itssigmonial & full(YESNO);
    i = i + 1;
end
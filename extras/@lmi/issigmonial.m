function itssigmonial=issigmonial(F)
%ISSIGMONIAL Check if all constraints are sigmonial       
%
% p = islinear(F)
%
% F : SET object
% p : boolean 0/1

monomtable = yalmip('monomtable');
xvars = [];

itssigmonial = 1;
i = 1;
F = flatten(F);
while itssigmonial & (i<=length(F.clauses))
    Fi = F.clauses{i};
    xvars = getvariables(Fi.data);
    monomtableX = monomtable(xvars,:);
    YESNO = any(find(any(0>monomtableX,2) | any(monomtableX-fix(monomtableX),2)));   
    itssigmonial = itssigmonial & full(YESNO);
    i = i + 1;
end
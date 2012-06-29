function [exponents,base]=getexponentbase(p,x)
%GETEXPONENTBASE Internal function used in SOS programs

% Author Johan Löfberg
% $Id: getexponentbase.m,v 1.9 2006-02-15 21:41:35 joloef Exp $

if isempty(p)
    exponents=[];
    base=[];
else
    p_vars = getvariables(p);
    x_vars = getvariables(x);

    base = getbase(p);

    monom_table = yalmip('monomtable');
    exponents = monom_table(p_vars,x_vars);
    if any(base(:,1))%base(1)~=0
        exponents = [spalloc(1,size(exponents,2),0);exponents];
    else
        base = base(:,2:end);
    end
end
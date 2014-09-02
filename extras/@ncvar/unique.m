function x = unique(x)
%UNIQUE (overloaded)

if min(x.dim(1),x.dim(2))~=1
    error('UNIQUE currently only supported for vectors')
end

base = x.basis;
vars = x.lmi_variables;
[new_base,i,j] = unique(base,'rows');

x.basis = new_base;
if x.dim(2)==1
    x.dim(1) = size(new_base,1);
else
    x.dim(2) = size(new_base,1);
end
% Reset info about conic terms
x.conicinfo = [0 0];


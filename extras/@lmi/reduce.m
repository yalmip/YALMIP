function F = projection(F,x,method)
% REDUCE  Removes redundant constraints using MPT
%
% Freduced  = reduce(F)
%
% F  :  Polytopic constraint object
%
% See also POLYTOPE, PROJECTION

f = [];
F = flatten(F);
for i = 1:length(F)
    if  F.clauses{i}.type==2
        fi =  F.clauses{i}.data;
        f = [f;fi(:)];
    else
        error('Only linear element-wise inequalities can be reduced')
    end
end

if ~islinear(F)
    error('Only linear element-wise inequalities can be reduced')
end

B = full(getbase(f));
P = polytope(-B(:,2:end),B(:,1));

x_vars = getvariables(F);
x = recover(x_vars);

H = get(P,'H');
K = get(P,'K');

F = (H*x <= K);

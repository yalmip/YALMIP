function [p,feasible,vol_reduction] = domain_reduction(p,upper,lower,lpsolver,xmin);
% This is just too expensive
LU = [p.lb p.ub];
t1 = p.binary_variables;
t2 = p.integer_variables;
p.binary_variables = [];
p.integer_variables = [];
if ~p.options.bmibnb.lpreduce | ((size(p.lpcuts,1)==0) & (any(p.lb(p.linears)<-1e8) & any(p.ub(p.linears)>1e8)))
    vol_reduction = 1;
    p.feasible = 1;

    p.lb(p.integer_variables) = ceil(p.lb(p.integer_variables));
    p.ub(p.integer_variables) = floor(p.ub(p.integer_variables));
    p.lb(p.binary_variables) = ceil(p.lb(p.binary_variables));
    p.ub(p.binary_variables) = floor(p.ub(p.binary_variables));

else
    [p,p.feasible] =  boxreduce(p,upper,lower,lpsolver,p.options,xmin);
end
p.binary_variables  = t1;
p.integer_variables = t2;
if ~isequal(LU,[p.lb p.ub])
    p.changedbounds = 1;
end
function [p,feasible,vol_reduction,seen_x] = domain_reduction(p,upper,lower,lpsolver,xmin);
LU = [p.lb p.ub];
seen_x = {};
if ~p.options.bmibnb.lpreduce | ((size(p.lpcuts,1)==0) & (any(p.lb(p.linears)<-1e8) & any(p.ub(p.linears)>1e8)))
    vol_reduction = 1;
    p.feasible = 1;
    p.lb(p.integer_variables) = ceil(p.lb(p.integer_variables)-1e-7);
    p.ub(p.integer_variables) = floor(p.ub(p.integer_variables)+1e-7);
    p.lb(p.binary_variables) = ceil(p.lb(p.binary_variables)-1e-7);
    p.ub(p.binary_variables) = floor(p.ub(p.binary_variables)+1e-7);
else
    [p,p.feasible,seen_x] =  boxreduce(p,upper,lower,lpsolver,p.options,xmin);
end
if ~isequal(LU,[p.lb p.ub])
    p.changedbounds = 1;
end
feasible = p.feasible;
vol_reduction = 0;
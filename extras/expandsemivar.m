function  [F_struc,K,binary_variables] = expandsemivar(F_struc,K,semicont_variables)

model.F_struc = F_struc;
model.K = K;
model.lb = -inf(size(model.F_struc,2)-1,1);
model.ub = -model.lb
model = presolve_bounds_from_modelbounds(model,1);

if any(isinf(model.lb(semicont_variables))) ||  any(isinf(model.ub(semicont_variables)))
    error('There are unbounded semi-continuous variables.');
end

% m new binaries required. 
m = length(semicont_variables);






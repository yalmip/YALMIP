function model = presolve_bounds_from_modelbounds(model);
if ~isempty(model.F_struc)
    [L,U] = findulb(model.F_struc,model.K);
    model.lb = max([model.lb L],[],2);
    model.ub = min([model.ub U],[],2);    
end



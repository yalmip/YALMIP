function x0 = create_trivial_initial(model)
x0 = (model.lb+model.ub)/2;
x0(isinf(model.ub)) = model.lb(isinf(model.ub))+1;
x0(isinf(model.lb)) = model.ub(isinf(model.lb))-1;
x0(isinf(x0)) = 0;
if any(model.variabletype == 4)
    problematic = find(any(model.monomtable(:,model.linearindicies) < 0 ,1));
    if ~isempty(problematic)
        problematic = problematic(find(x0(problematic)==0));
        Oneisfeas = problematic(find(model.ub(problematic) > 1));
        x0(Oneisfeas) = 1;
    end
end
x0(find(model.lb==model.ub)) = model.lb(find(model.lb==model.ub));
end
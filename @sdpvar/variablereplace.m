function y = variablereplace(y,oldVar,newVar)

var = y.lmi_variables;

[~,pos] = ismember(oldVar,var);
if any(pos)
    index = find(pos);
    pos = pos(index);
    y.lmi_variables(pos) = newVar(index);
    
    [var,pos] = sort(y.lmi_variables);
    y.lmi_variables = var;
    y.basis = [y.basis(:,1) y.basis(:,1+pos)];
end

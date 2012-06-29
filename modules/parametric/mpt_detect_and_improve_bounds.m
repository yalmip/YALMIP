function [lb,ub,redundant] = detect_and_improve_bounds(Matrices,lb,ub,binary_var_index,options);

A = [ Matrices.G -Matrices.E];
b = [ Matrices.W ];

[global_lower,global_upper,bound_rows] = find_variable_bounds(A,b,[Matrices.Aeq Matrices.Beq],Matrices.beq);
global_lower(binary_var_index) = max(global_lower(binary_var_index),0);
global_upper(binary_var_index) = min(global_upper(binary_var_index),1);

if ~isempty(lb)
    global_lower = max([global_lower lb],[],2);
end
if ~isempty(ub)
    global_upper = min([global_upper ub],[],2);
end

[lb,ub,redundant,psstruct,infeasible] = tightenbounds(A,b,global_lower,global_upper,[],binary_var_index);

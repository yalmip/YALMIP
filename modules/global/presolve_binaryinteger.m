function p = presolve_binaryinteger(p)
% Promote integers to binary if possible
if ~isempty(p.integer_variables)
    binary_integer = find(p.ub(p.integer_variables)<=1 & p.lb(p.integer_variables)>=0);
    if ~isempty(binary_integer)
        p.binary_variables = union(p.binary_variables,p.integer_variables(binary_integer));
        p.isinteger(p.integer_variables(binary_integer)) = 0;
        p.integer_variables(binary_integer) = [];
        p.isbinary(p.binary_variables) = 1;
    end
end
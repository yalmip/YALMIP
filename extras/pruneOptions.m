function options = pruneOptions(options)

if isequal(options.mosek,options.default.mosek)
    options.mosek = [];
end
if isequal(options.cplex,options.default.cplex)
    options.mosek = [];
end
if isequal(options.gurobi,options.default.gurobi)
    options.mosek = [];
end
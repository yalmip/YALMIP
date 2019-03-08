function  p = adjustMaxTime(p,maxtime,elapsed)

% Allow the solver at-least 1 second
remaining = maxtime-elapsed;
budgetForSolver = max(1,ceil(remaining));

% TODO: Support more solvers
p.options.cplex.timelimit = budgetForSolver;
p.options.mosek.MSK_DPAR_MIO_MAX_TIME = budgetForSolver;
p.options.gurobi.TimeLimit = budgetForSolver;

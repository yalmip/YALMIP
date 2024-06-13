function [solver,diagnostic] = setupBNB(solver,ProblemClass,options,solvers,socp_are_really_qc,F,h,logdetStruct,parametric,evaluation_based,F_vars,exponential_cone,allsolvers)

diagnostic = [];
options.forceglobal = 0;
temp_options = options;
temp_options.solver = options.bnb.solver;
tempProblemClass = ProblemClass;
tempProblemClass.constraint.binary  = 0;
tempProblemClass.constraint.integer = 0;
tempProblemClass.constraint.semicont = 0;
localsolver = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
if isempty(localsolver) | strcmpi(localsolver.tag,'bnb')
    if isempty(temp_options.bnb.solver)
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-2,'YALMIP');
        diagnostic.problem = -2;
    else
        % User has specified a lower-bound solver, but we failed to use
        % this. Could be that it doesn't exist, or that it is not
        % applicable
        for i = 1:length(solvers)
            if strcmpi(solvers(i).tag,temp_options.bnb.solver)
                % Solver exist, hence it is not applicable
                diagnostic.solvertime = 0;
                diagnostic.info = yalmiperror(-4,temp_options.bnb.solver);
                diagnostic.problem = -4;
                return
            end
        end
        % Solver was not found in list of available solvers, hance user
        % has specified a solver which doesn't exist
        diagnostic.solvertime = 0;
        diagnostic.info = yalmiperror(-3,temp_options.bnb.solver);
        diagnostic.problem = -3;
        return
    end
    return
elseif strcmpi(localsolver.tag,'bmibnb')
    [localsolver,diagnostics] = setupBMIBNB(localsolver,tempProblemClass,options,solvers,socp_are_really_qc,F,h,logdetStruct,parametric,evaluation_based,F_vars,[],allsolvers);
end
solver.lower = localsolver;

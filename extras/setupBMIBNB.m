function [solver,diagnostic] = setupBMIBNB(solver,ProblemClass,options,solvers,socp_are_really_qc,F,h,logdetStruct,parametric,evaluation_based,F_vars,exponential_cone,allsolvers)

diagnostic = [];

% Relax problem for lower solver
tempProblemClass = ProblemClass;

sdp = tempProblemClass.constraint.inequalities.semidefinite;
tempProblemClass.constraint.inequalities.semidefinite.linear = sdp.linear | sdp.quadratic | sdp.polynomial;
tempProblemClass.constraint.inequalities.semidefinite.quadratic = 0;
tempProblemClass.constraint.inequalities.semidefinite.polynomial = 0;
tempProblemClass.constraint.inequalities.rank = 0;

lp = tempProblemClass.constraint.inequalities.elementwise;
tempProblemClass.constraint.inequalities.elementwise.linear = lp.linear | lp.quadratic.convex | lp.quadratic.nonconvex | sdp.polynomial;
tempProblemClass.constraint.inequalities.elementwise.quadratic.convex = 0;
tempProblemClass.constraint.inequalities.elementwise.quadratic.nonconvex = 0;
tempProblemClass.constraint.inequalities.elementwise.polynomial = 0;
tempProblemClass.constraint.inequalities.elementwise.sigmonial = 0;

equ = tempProblemClass.constraint.equalities;
tempProblemClass.constraint.equalities.linear = equ.linear | equ.quadratic | equ.polynomial;
tempProblemClass.constraint.equalities.quadratic = 0;
tempProblemClass.constraint.equalities.polynomial = 0;
tempProblemClass.constraint.equalities.sigmonial = 0;

tempProblemClass.objective.quadratic.nonconvex = 0;
tempProblemClass.objective.polynomial = 0;
tempProblemClass.objective.sigmonial = 0;

tempProblemClass.constraint.inequalities.rank  = 0;
tempProblemClass.evaluation  = 0;
tempProblemClass.exponentialcone  = 0;

temp_options = options;
temp_options.solver = options.bmibnb.lowersolver;

% If the problem actually is quadratic, try to get a convex problem
% this will typically allow us to solver better lower bounding problems
% (we don't have to linearize the cost)
[lowersolver,problem] = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
if isempty(lowersolver) || strcmpi(lowersolver.tag,'bmibnb') || strcmpi(lowersolver.tag,'bnb')
    % No, probably non-convex cost. Pick a linear solver instead and go
    % for lower bound based on a complete "linearization"
    tempProblemClass.objective.quadratic.convex = 0;
    [lowersolver,problem] = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
end

if isempty(lowersolver) || strcmpi(lowersolver.tag,'bmibnb') || strcmpi(lowersolver.tag,'bnb')
    tempbinary = tempProblemClass.constraint.binary;
    tempinteger = tempProblemClass.constraint.integer;
    tempsemicont = tempProblemClass.constraint.semicont;
    tempProblemClass.constraint.binary = 0;
    tempProblemClass.constraint.integer = 0;
    tempProblemClass.constraint.semicont = 0;
    [lowersolver,problem] = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    tempProblemClass.constraint.binary = tempbinary;
    tempProblemClass.constraint.integer = tempinteger;
    tempProblemClass.constraint.semicont = tempsemicont;
end

if isempty(lowersolver) || strcmpi(lowersolver.tag,'bmibnb') || strcmpi(lowersolver.tag,'bnb')
    diagnostic.solvertime = 0;
    diagnostic.info = yalmiperror(-2,'YALMIP');
    diagnostic.problem = -2;
    return
end
solver.lowercall = lowersolver.call;
solver.lowersolver = lowersolver;

temp_options = options;
temp_options.solver = options.bmibnb.uppersolver;

if any(getcutflag(F))
    % Could be that the model involves, e.g., semidefinite cuts, which
    % shouldn't be sent to the upper bound solver
    Ftemp = F;
    Ftemp(find(getcutflag(F)))=[];
    [temp_ProblemClass,aux1,aux2,aux3,aux4,aux5,aux6] = categorizeproblem(Ftemp,logdetStruct,h,options.relax,parametric,evaluation_based,F_vars,exponential_cone);
    temp_ProblemClass.gppossible = ProblemClass.gppossible;
else
    temp_ProblemClass = ProblemClass;
end

temp_ProblemClass.constraint.binary = 0;
temp_ProblemClass.constraint.integer = 0;
[uppersolver,problem] = selectsolver(temp_options,temp_ProblemClass,solvers,socp_are_really_qc,allsolvers);
if ~isempty(uppersolver) && strcmpi(uppersolver.tag,'bnb')
    temp_options.solver = 'none';
    [uppersolver,problem] = selectsolver(temp_options,temp_ProblemClass,solvers,socp_are_really_qc,allsolvers);
end
if isempty(uppersolver) || strcmpi(uppersolver.tag,'bmibnb')
    diagnostic.solvertime = 0;
    diagnostic.info = yalmiperror(-2,'YALMIP');
    diagnostic.problem = -2;
    return
end
if strcmpi(uppersolver.version,'geometric') &&  strcmpi(uppersolver.tag,'fmincon')
    uppersolver.version = 'standard';
    uppersolver.call = 'callfmincon';
end
if strcmpi(uppersolver.version,'geometric') &&  strcmpi(uppersolver.tag,'ipopt')
    uppersolver.version = 'standard';
    uppersolver.call = 'callipoptmex';
end
if strcmpi(uppersolver.version,'geometric') &&  strcmpi(uppersolver.tag,'snopt')
    uppersolver.version = 'standard';
    uppersolver.call = 'callsnopt';
end
if strcmpi(uppersolver.version,'geometric') &&  strcmpi(uppersolver.tag,'pennon')
    uppersolver.version = 'standard';
    uppersolver.call = 'callpennonm';
end

solver.uppercall = uppersolver.call;
solver.uppersolver = uppersolver;

temp_options = options;
temp_options.solver = options.bmibnb.lpsolver;
tempProblemClass.constraint.inequalities.semidefinite.linear = 0;
tempProblemClass.constraint.inequalities.semidefinite.quadratic = 0;
tempProblemClass.constraint.inequalities.semidefinite.polynomial = 0;
tempProblemClass.constraint.inequalities.secondordercone.linear = 0;
tempProblemClass.constraint.inequalities.secondordercone.nonlinear = 0;
tempProblemClass.objective.quadratic.convex = 0;
tempProblemClass.objective.quadratic.nonconvex = 0;
tempProblemClass.objective.quadratic.nonconvex = 0;
tempProblemClass.objective.polynomial = 0;
tempProblemClass.objective.sigmonial = 0;

[lpsolver,problem] = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);

if isempty(lowersolver) || strcmpi(lowersolver.tag,'bmibnb')
    tempbinary = tempProblemClass.constraint.binary;
    tempProblemClass.constraint.binary = 0;
    [lpsolver,problem] = selectsolver(temp_options,tempProblemClass,solvers,socp_are_really_qc,allsolvers);
    tempProblemClass.constraint.binary = tempbinary;
end

if isempty(lpsolver) || strcmpi(lpsolver.tag,'bmibnb')
    diagnostic.solvertime = 0;
    diagnostic.info = yalmiperror(-2,'YALMIP');
    diagnostic.problem = -2;
    return
end
solver.lpsolver = lpsolver;
solver.lpcall = lpsolver.call;
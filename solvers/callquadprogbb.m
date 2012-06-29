function output = callquadprog(interfacedata)

% Author Johan Löfberg
% $Id: callquadprog.m,v 1.17 2007-08-02 11:39:36 joloef Exp $

options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
solveroutput = callsolver(model,options);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end
solution = quadprogsol2yalmipsol(solveroutput,model);

% Save all data sent to solver?
if ~options.savesolverinput
    model = [];
end

% Save all data from the solver?
if ~options.savesolveroutput
    solveroutput = [];
end

% Standard interface
output.Primal      = solution.x(:);
output.Dual        = solution.D_struc;
output.Slack       = [];
output.problem     = solution.problem;
output.infostr     = yalmiperror(solution.problem,interfacedata.solver.tag);
output.solvertime  = solvertime;
if ~options.savesolverinput
    output.solverinput = [];
else
    output.solverinput = model;
end
if ~options.savesolveroutput
    output.solveroutput = [];
else
    output.solveroutput = solveroutput;
end

function solveroutput = callsolver(model,options)
x = [];
fmin = [];
flag = [];
output = [];
lambda = [];
options.quadprogbb.constant = model.f;
if options.verbose
    options.quadprogbb.verbosity = options.verbose - 1;
    [x,fval,time,stat] = quadprogbb(model.Q,model.c,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,options.quadprogbb); 
else
    options.quadprogbb.verbosity = 0;
    evalc('[x,fval,time,stat] = quadprogbb(model.Q,model.c,model.A,model.b,model.Aeq,model.beq,model.lb,model.ub,options.quadprogbb);');
end
solveroutput.x = x;
solveroutput.fval = fval;
solveroutput.time = time;
solveroutput.stat = stat;

function solution = quadprogsol2yalmipsol(solveroutput,model)

solution.x = solveroutput.x(:);
solution.D_struc = [];
switch solveroutput.stat.status
    case 'opt_soln'
        solution.problem = 0;
    case 'inf_or_unb'
        solution.problem = 12;
    case 'num_issues'
        solution.problem = 4;
    case 'time_limit'
        solution.problem = 3;
    otherwise
        solution.problem = -1;
end
function output = callqpip(interfacedata)

% Author Johan Löfberg
% $Id: callqpip.m,v 1.4 2007-05-22 13:40:36 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if options.savedebug
    model.ops = options.qpip;
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = clock;
[x,flag,lm] = qpip(model.Q, model.c, model.A, model.b, model.Aeq, model.beq, model.lb, model.ub, options.verbose,options.qpip.mu,options.qpip.method);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Internal format for duals
if ~isempty(lm)
    D_struc = [lm.equality;lm.inequality];
else
    D_struc = [];
end

% Check, currently not exhaustive...
switch flag
    case 0
        problem = 0;
    case 2
        problem = 1;
    otherwise
        problem = 1;
end

infostr = yalmiperror(problem,'QPIP');

% Save all data sent to solver?
if options.savesolverinput
    solverinput = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.flag = flag;
    solveroutput.lm = lm;
else
    solveroutput = [];
end

% Standard interface
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;
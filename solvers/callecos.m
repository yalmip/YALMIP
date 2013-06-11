function output = callecos(yalmipmodel)

% Retrieve needed data
options = yalmipmodel.options;

model = yalmip2ecos(yalmipmodel);

if options.savedebug
    save ecosdebug model
end

if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end
solvertime = clock;
if isempty(model.A)
    if options.verbose
        [x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims);
    else
        evalc('[x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims);');
    end
else
    if options.verbose
        [x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims,model.A,model.b);
    else
        evalc('[x,y,info,s,z] = ecos(model.c,model.G,model.h,model.dims,model.A,model.b);');
    end
end
if yalmipmodel.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Internal format
Primal = x;
Dual   = [y;z];

switch info.exitflag
    case 0
        problem = 0;
    case 1
        problem = 1;
    case 2
        problem = 2;
    case -1
        problem = 3;
    case {-2,-3}
        problem = 4;
    case -7
        problem = 9;
    otherwise
        problem = 9;
end

infostr = yalmiperror(problem,yalmipmodel.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.y = y;
    solveroutput.info = info;
    solveroutput.s = s;
    solveroutput.z = z;
else
    solveroutput = [];
end

% Standard interface
output.Primal      = Primal;
output.Dual        = Dual;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;
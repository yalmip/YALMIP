function output = callsedumi(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
K       = model.K;
ub      = model.ub;
lb      = model.lb;

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(ub)
    [F_struc,K] = addbounds(F_struc,K,ub,lb);
end

data.A = -F_struc(:,2:end);
data.b = full(F_struc(:,1));
data.c = c;
cones = K;

if options.savedebug
    save scsdebug data cones
end

if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end

solvertime = clock; 
problem = 0;  
params = options.scs;
params.verbose = options.verbose;
switch  model.solver.tag
    case 'scs-direct'
         [x_s,y_s,s,info] = scs_direct(data,cones,params);
    otherwise
        [x_s,y_s,s,info] = scs_indirect(data,cones,params);
end

% solvertime = cputime - solvertime;%etime(clock,solvertime
if model.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Internal format
Primal = x_s; 
Dual   = y_s;

switch info.status
    case 'Solved'
        problem = 0;
    case 'Infeasible'
        problem = 1;
     case 'Unbounded'
        problem = 2;    
    otherwise
        status = 9;
end

infostr = yalmiperror(problem,model.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.data = data
    solverinput.cones = cones;
    solverinput.param = param;  
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.x = x_s;
    solveroutput.y = y_s;
    solveroutput.s = s;
    solveroutput.info = info;
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
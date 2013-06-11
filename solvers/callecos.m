function output = callecos(model)

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

if K.f > 0
    b = -full(F_struc(1:K.f,1));
    A = F_struc(1:K.f,2:end);
    h = full(F_struc(1+K.f:end,1));
    G = -F_struc(1+K.f:end,2:end);
else
    A = [];
    b = [];
    h = full(F_struc(1:end,1));
    G = -F_struc(1:end,2:end);
end

dims.l = K.l;
if K.q(1)==0
    dims.q = [];
else
    dims.q = K.q;
end

if options.savedebug
    save ecosdebug  A b G h c dims
end

if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end
solvertime = clock;
if isempty(A)
    if options.verbose
        [x,y,info,s,z] = ecos(c,G,h,dims);
    else
        evalc('[x,y,info,s,z] = ecos(c,G,h,dims);');
    end
else
    if options.verbose
        [x,y,info,s,z] = ecos(c,G,h,dims,A,b);
    else
        evalc('[x,y,info,s,z] = ecos(c,G,h,dims,A,b);');
    end
end
if model.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

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
    otherwise
        problem = 9;
end

infostr = yalmiperror(problem,model.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.G = G;
    solverinput.h = h;
    solverinput.c = c;
    solverinput.dims = dims;
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
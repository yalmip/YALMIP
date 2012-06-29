function output = calloptiqsopt(interfacedata)

% Author Johan Löfberg 
% $Id: callopticlp.m,v 1.6 2005-05-07 13:53:20 joloef Exp $

% Standard input interface
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

showprogress('Calling QSOPT',options.showprogress);

if isempty(F_struc)  
    A = [];
    b = [];
else
    A =[-F_struc(K.f+1:end,2:end);-F_struc(1:K.f,2:end)];
    b =[F_struc(K.f+1:end,1);F_struc(1:K.f,1)];    
end

opts.nin = length(b)-K.f;
%opts.tolfun = options.clp.primaltolerance;
opts.maxiter = options.clp.maxnumiterations;
opts.maxtime = options.clp.maxnumseconds;
opts.display = options.verbose;

if options.savedebug    
    save qsoptdebug c A b  lb ub opts
end

solvertime = clock; 
[x,fval,exitflag] = qsopt(full(c), A, full(b), full(lb), full(ub),opts);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% No duals
D_struc = [];

switch exitflag
    case 1
        problem = 0;
    case 0
        problem = 3;
    case -1
        problem = 1;
    case -2
        problem = 11;   
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'QSOPT');       

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.f = c;
    solverinput.lb = lb;
    solverinput.ub = ub;    
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
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
function output = calloptiqsopt(interfacedata)

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

opts = options.qsopt;
opts.nin = length(b)-K.f;
opts.display = options.verbose;

if options.savedebug    
    save qsoptdebug c A b  lb ub opts
end

solvertime = tic;
[x,fval,exitflag] = qsopt(full(c), A, full(b), full(lb), full(ub),opts);
solvertime = toc(solvertime);

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
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
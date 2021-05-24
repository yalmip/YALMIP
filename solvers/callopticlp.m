function output = callopticlp(interfacedata)

% Standard input interface
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
Q       = interfacedata.Q;

showprogress('Calling CLP',options.showprogress);

if isempty(F_struc)  
    A = sparse([]);
    b = [];
else
    A = -F_struc(1:K.f+K.l,2:end);
    b = F_struc(1:K.f+K.l,1);    
end

opts = options.clp;

if length(b)>0
    rl = repmat(-inf,length(b),1);
    rl(1:K.f) = b(1:K.f);
else
    rl = [];
end
ru = full(b);
lb = full(lb);
ub = full(ub);

H = 2*sparse(tril(Q));
if options.savedebug
    save clpdebug c A b  lb ub opts H
end

solvertime = tic;
[x,fval,exitflag,iter,lambda] = clp(H,full(c), A, rl, ru, lb, ub,opts);
solvertime = toc(solvertime);

% No duals
D_struc = -lambda.dual_row;

switch exitflag
    case 0
        problem = 0;
    case 1
        problem = 1;
    case 2
        problem = 2;
    case 3
        problem = 3;
    case 5
        problem = 16;
    case 4
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
    solveroutput.iter = iter;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x(:),D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
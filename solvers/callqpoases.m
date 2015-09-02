function output = callqpoases(interfacedata)

options = interfacedata.options;
model = yalmip2quadprog(interfacedata);

if options.savedebug
    save debugfile model
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
solvertime = tic;
solution = callsolver(model,options);
solvertime = toc(solvertime);

switch solution.exitflag
    case 0
        problem = 0;
    case 1
        problem = 3;
    case -2
        problem = 1;
    case -3
        problem = 2;
    case -1
        problem = 11;
    otherwise
        problem = -1;
end

% Standard interface
Primal      = solution.x(:);
Dual        = solution.lambda(:);
infostr     = yalmiperror(problem,interfacedata.solver.tag);
if ~options.savesolverinput
    solverinput = [];
else
    solverinput = model;
end
if ~options.savesolveroutput
    solveroutput = [];
else
    solveroutput = solution;
end
% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,infostr,solverinput,solveroutput,solvertime);


function solveroutput = callsolver(model,options)
x = [];
fval = [];
exitflag = [];
iter = [];
lambda = [];
lbA = full([model.beq;-inf(length(model.b),1)]);
ubA = full([model.beq;model.b]);
A = [model.Aeq;model.A];
if nnz(model.Q) == 0
    options.qpoases.enableRegularisation=1;
end
options.qpoases.printLevel = options.verbose+1;
[x,fval,exitflag,iter,lambda] = qpOASES(model.Q, model.c, A, model.lb,model.ub,lbA,ubA,options.qpoases,qpOASES_auxInput());
solveroutput.x = x;
solveroutput.fval = fval;
solveroutput.exitflag = exitflag;
solveroutput.iter = iter;
if ~isempty(lambda)
    solveroutput.lambda = -lambda(1+length(model.lb):end);
else
    solveroutput.lambda=[];
end

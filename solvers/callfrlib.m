function output = callfrlib(model)

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(model.ub)
    [model.F_struc,model.K] = addStructureBounds(model.F_struc,model.K,model.ub,model.lb);
end

options = model.options;

% SeDuMi format
A = -model.F_struc(:,2:end)';
b = -model.c;
C = model.F_struc(:,1);
K = model.K;

if (strcmp(options.frlib.reduce,'') || strcmp(options.frlib.reduce,'auto'))
    if model.dualized
        options.frlib.reduce = 'primal';
    else
        options.frlib.reduce = 'dual';
    end
end

if options.savedebug
    save frlibdebug A b C K
end

% *********************************************
% Call frlib and solver
% *********************************************
if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end

% Reduce
prg = frlibPrg(A,b,C,K);

switch options.frlib.reduce
    case 'primal'
        prgR = prg.ReducePrimal(options.frlib.approximation,options.frlib);
    case 'dual'
        prgR = prg.ReduceDual(options.frlib.approximation,options.frlib);
    otherwise
       error('Wrong option ''reduce'' for frlib.');
end

if options.verbose
    PrintStats(prgR);
end

% Back to YALMIP internal format to enable arbitrary solver
model.F_struc = [prgR.c' -prgR.A'];
model.c = -prgR.b;
model.K = prgR.K;
if any(model.K.s) && any(model.K.s == 0)
    model.K.s(model.K.s==0)=[];
end
if any(model.K.q) && any(model.K.q == 0)
    model.K.q(model.K.q==0)=[];
end

% Call SDP solver
solvertime = tic;
output = feval(model.solver.solver.call,model);
solvertime = toc(solvertime);
 
[x,y,dual_recov_success] = prgR.Recover(output.Dual,output.Primal);
output.Dual = x;
output.Primal = y;

output.infostr = yalmiperror(output.problem,[model.solver.tag ' -> ' model.solver.solver.tag]);

if options.savesolverinput   
    output.solverinput.A = A';
    output.solverinput.b = b;
    output.solverinput.c = C;
    output.solverinput.K = K;
end
if options.savesolveroutput  
    output.solveroutput.reducedmodel = prgR;
    output.solveroutput.dual_recov_success = dual_recov_success;
end

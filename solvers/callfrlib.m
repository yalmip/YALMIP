function output = callfrlib(model)

% *********************************************
% Bounded variables converted to constraints
% N.B. Only happens when caller is BNB
% *********************************************
if ~isempty(model.ub)
    [model.F_struc,model.K] = addbounds(model.F_struc,model.K,model.ub,model.lb);
end

options = model.options;

% SeDuMi format
A = -model.F_struc(:,2:end)';
b = -model.c;
C = model.F_struc(:,1);
K = model.K;

if options.savedebug
    save frlibdebug A b C K
end

% *********************************************
% Call frlib and solver
% *********************************************
if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end
solvertime = clock; 
% Reduce
prg = frlibPrg(A,b,C,K);
switch options.frlib.reduce
    case 'primal'
        prgR = prg.ReducePrimal(options.frlib.approximation);
    case 'dual'
        prgR = prg.ReduceDual(options.frlib.approximation);
    otherwise
end
if options.verbose
    PrintStats(prgR);
end
% Use frlibs clean to remove silly stuff
[A,b,T] = CleanLinear(prgR.A,prgR.b,options.frlib.useQR);

% Back to YALMIP internal format to enable arbitrary solver
model.F_struc = [prgR.c' -A'];
model.c = -b;
model.K = prgR.K;
if any(model.K.s) && any(model.K.s == 0)
    model.K.s(model.K.s==0)=[];
end
if any(model.K.q) && any(model.K.q == 0)
    model.K.q(model.K.q==0)=[];
end

% Call SDP solver
output = feval(model.solver.solver.call,model);

[x,y,dual_recov_success] = prgR.Recover(output.Dual,T*output.Primal);
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
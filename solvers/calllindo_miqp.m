function output = calllindo_miqp(interfacedata)

global MY_LICENSE_FILE
lindo

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;
monomtable = interfacedata.monomtable;

lindo;

if ~isempty(F_struc)
    A = -F_struc(:,2:end);
    b = F_struc(:,1);   
    csense = [repmat('E',1,K.f) repmat('L',1,K.l)];
else
    A = [];
    b = [];
    csense = [];
    
    A = ones(1,length(c));
    b = 1e6;
    csense = 'L';
end

% Specifying variable types...
vtype = repmat('C',1,length(c));
vtype(interfacedata.integer_variables) = 'I';

% Specifying the quadratic portion of the problem data 
[QCvar1,QCvar2,QCcoef] = find(triu(Q));
QCrows = zeros(1,length(QCvar1));
QCvar1 = QCvar1 - 1;
QCvar2 = QCvar2 - 1;

% Solve the problem using the generic QP/LP/MIP/MIQP solver (lmsolvemp.m)
objsen=LS_MIN;
solver=eval(options.lindo.LS_METHOD);

solvertime = tic;
if nnz(Q)>0
    [x,D_struc,s,dj,pobj,solstat,nErr]  = LMsolvem(A,full(b),full(c),csense,lb,ub,vtype,QCrows-1,QCvar1,QCvar2,2*QCcoef,objsen,solver,options.verbose);
else
    [x,D_struc,s,dj,pobj,solstat,nErr]  = LMsolvem(A,full(b),full(c),csense,lb,ub,vtype,[],[],[],[],objsen,solver,options.verbose);
end
solvertime = toc(solvertime);
    
switch solstat
    case {LS_STATUS_OPTIMAL,LS_STATUS_BASIC_OPTIMAL,7,8}
        problem = 0;
    case {LS_STATUS_INFEASIBLE}
        problem = 1;
    case {LS_STATUS_UNBOUNDED}
        problem = 2;
    otherwise
        problem = 11;
end
infostr = yalmiperror(problem,'LINDO-QP');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.A = A;
    solverinput.b = b;
    solverinput.c = c;
    solverinput.csense = csense;
    solverinput.vtype = vtype;    
    solverinput.lb = lb;    
    solverinput.vtype = vtype;    
    solverinput.objsen = objsen;   
    solverinput.solver = solver;   
    solverinput.options = options.fmincon;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.w = x;
    solveroutput.y = D_struc;
    solveroutput.s = s;
    solveroutput.pobj=pobj;
    solveroutput.solstat=solstat;
    solveroutput.nErr=nErr;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(x,D_struc,[],problem,'LINDO',solverinput,solveroutput,solvertime);
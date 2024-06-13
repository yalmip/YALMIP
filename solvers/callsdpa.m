function output = callsdpa(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

% Bounded variables converted to constraints
if ~isempty(ub)
    [F_struc,K] = addStructureBounds(F_struc,K,ub,lb);
end

% Convert from internal (sedumi) format
[mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(F_struc,c,K);

if options.verbose==0
    options.sdpa.print = 'no';
else
    options.sdpa.print = 'display';
end

if options.savedebug
    ops = options.sdpa;
    save sdpadebug mDIM nBLOCK bLOCKsTRUCT c F ops
end

if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end

solvertime = tic;
[objVal,x,X,Y,INFO]=sdpam(mDIM,nBLOCK,bLOCKsTRUCT,c,F,[],[],[],options.sdpa);
solvertime = toc(solvertime);

% Create variables in YALMIP internal format
Primal = x;

Dual = [];
for i = 1:length(Y)
    Dual = [Dual;Y{i}(:)];
end

Slack = [];
if options.saveduals
    for i = 1:length(X)
        Slack = [Slack;X{i}(:)];
    end
end

switch (INFO.phasevalue)
    case 'pdOPT'
        problem = 0;
    case {'noINFO','pFEAS','dFEAS'}
        problem = 3;
    case {'pdFEAS'}
        problem = 4;
    case 'pFEAS_dINF'
        problem = 2;
    case 'pINF_dFEAS'
        problem = 1;
    case 'pUNBD'
        problem = 2;
    case 'dUNBD'
        problem = 1;
    case 'pdINF'
        problem = 12;
    otherwise
        problem = -1;
end

if options.savesolveroutput
    solveroutput.objVal = objVal;
    solveroutput.x = x;
    solveroutput.X = X;
    solveroutput.Y = Y;
    solveroutput.INFO = INFO;
else
    solveroutput = [];
end

if options.savesolverinput
    solverinput.mDIM = mDIM;
    solverinput.nBLOCK=nBLOCK;
    solverinput.bLOCKsTRUCT=bLOCKsTRUCT;
    solverinput.c=c;
    solverinput.F=F;
else
    solverinput = [];
end

% Standard interface 
output = createOutputStructure(Primal,Dual,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
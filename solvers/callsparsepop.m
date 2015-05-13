function output = callsparsepop(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
f       = model.f;
K       = model.K;
ub      = model.ub;
lb      = model.lb;
monoms  = model.monomtable;

nvars = nnz(model.variabletype == 0);
monoms = monoms(:,model.variabletype == 0);

% Sparsepop exploit bounds
[lb,ub,cand_rows_eq,cand_rows_lp] = findlinearulb(F_struc,K,lb,ub,find(model.variabletype == 0));
lb = lb(find(model.variabletype==0));lb = lb(:)';
ub = ub(find(model.variabletype==0));ub = ub(:)';
F_struc(K.f+cand_rows_lp,:)=[];
F_struc(cand_rows_eq,:)=[];
K.l = K.l-length(cand_rows_lp);
K.f = K.f-length(cand_rows_eq);

obj.typeCone = 1;
obj.sizeCone = 1;
constraint = [];
if nnz(c)==0
    obj.degree = 0;
    obj.noTerms = 1;
    obj.dimVar = nvars;
    obj.coef     = f;
    obj.supports = [zeros(1,size(monoms,2))];
else
    obj.degree   = full(max(sum(monoms(find(c),:),2)));
    obj.noTerms  = nnz(c)+nnz(f);
    obj.dimVar   = nvars;
    obj.coef     = [f(find(f));c(find(c))];
    obj.supports = monoms(find(c),:);
    if find(f)
        obj.supports = [zeros(1,size(monoms,2)); obj.supports];
    end
end

for i = 1:K.f
    f0 = F_struc(i,1);
    f = F_struc(i,2:end);
    constraint{i}.typeCone = -1;
    constraint{i}.degree = full(max(sum(monoms(find(f),:),2)));
    constraint{i}.noTerms = nnz(f)+nnz(f0);
    constraint{i}.dimVar = nvars;
    constraint{i}.coef   = [f0(find(f0));f(find(f))'];
    constraint{i}.supports = monoms(find(f),:);
    if find(f0)
        constraint{i}.supports = [zeros(1,size(monoms,2)); constraint{i}.supports];
    end
end

for i = 1:K.l
    f0 = F_struc(i+K.f,1);
    f = F_struc(i+K.f,2:end);
    constraint{i+K.f}.typeCone = 1;
    constraint{i+K.f}.degree = full(max(sum(monoms(find(f),:),2)));
    constraint{i+K.f}.noTerms = nnz(f)+nnz(f0);
    constraint{i+K.f}.dimVar = nvars;
    constraint{i+K.f}.coef   = [f0(find(f0));f(find(f))'];
    constraint{i+K.f}.supports = monoms(find(f),:);
    if find(f0)
        constraint{i+K.f}.supports = [zeros(1,size(monoms,2)); constraint{i+K.f}.supports];
    end
end

paramin = options.sparsepop;

switch options.verbose
    case 0
        paramin.printLevel = [0 0];
    case 1
        paramin.printLevel = [1 0];
    case 2
        paramin.printLevel = [2 2];
end
if options.savedebug
    save sparsepopdebug obj constraint paramin
end

% *********************************************
% Call sparsePOP
% *********************************************
if options.showprogress;showprogress(['Calling ' model.solver.tag],options.showprogress);end
problem = 0;  
lb(isinf(lb)) = -1.0e10;
ub(isinf(ub)) = 1.0e10;

solvertime = tic;
if options.verbose==0
    evalc('[param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo] = sparsePOP(obj,constraint,lb,ub,paramin);');
else
    [param,SDPobjValue,POP,cpuTime,SDPsolverInfo,SDPinfo] = sparsePOP(obj,constraint,lb,ub,paramin);
end
solvertime = toc(solvertime);

if ~isempty(POP.xVect)
    Primal = zeros(length(c),1);
    Primal(model.variabletype==0) = POP.xVect;
else
    Primal = [];
end

infostr = yalmiperror(problem,model.solver.tag);

% Save ALL data sent to solver
if options.savesolverinput
    solverinput.obj = obj;
    solverinput.constraint = constraint;
    solverinput.lb = lb;
    solverinput.ub = ub;    
else
    solverinput = [];
end

% Save ALL data from the solution?
if options.savesolveroutput
    solveroutput.param =param;
    solveroutput.SDPobjValue = SDPobjValue;
    solveroutput.POP = POP;
    solveroutput.cpuTime = cpuTime;
    solveroutput.SDPsolverInfo = SDPsolverInfo;
    solveroutput.SDPinfo = SDPinfo;
else
    solveroutput = [];
end

% Standard interface
output = createOutputStructure(Primal,[],[],problem,infostr,solverinput,solveroutput,solvertime);
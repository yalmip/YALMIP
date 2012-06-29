function output = calllindo_nlp(interfacedata)

global MY_LICENSE_FILE

persistent iEnv

% Instead of calling lindo, we define the parameters we need. This is
% needed to speed up repeated calls
%lindo
LSERR_NO_ERROR                                               = 0000;
LS_IPARAM_NLP_PRINTLEVEL                                     = 203;
LS_IPARAM_NLP_SOLVER                                         = 201;
LS_IPARAM_NLP_MAXLOCALSEARCH                                 = 221;
LS_STATUS_OPTIMAL                                            = 1;
LS_STATUS_BASIC_OPTIMAL                                      = 2;
LS_STATUS_INFEASIBLE                                         = 3;
LS_STATUS_UNBOUNDED                                          = 4;
LS_STATUS_FEASIBLE                                           = 5;
LS_STATUS_INFORUNB                                           = 6;
LS_STATUS_NEAR_OPTIMAL                                       = 7;
LS_STATUS_LOCAL_OPTIMAL                                      = 8;
LS_STATUS_LOCAL_INFEASIBLE                                   = 9;
LS_STATUS_CUTOFF                                             = 10;
LS_STATUS_NUMERICAL_ERROR                                    = 11;
LS_STATUS_UNKNOWN                                            = 12;
LS_STATUS_UNLOADED                                           = 13;
LS_STATUS_LOADED                                             = 14;
LS_METHOD_FREE                                               = 0;
LS_METHOD_PSIMPLEX                                           = 1;
LS_METHOD_DSIMPLEX                                           = 2;
LS_METHOD_BARRIER                                            = 3;
LS_METHOD_NLP                                                = 4;
LS_NMETHOD_FREE                                              = 4;
LS_NMETHOD_CONOPT                                            = 7;
LS_NMETHOD_MSW_GRG                                           = 9;

if isempty(iEnv)
    % This call is mighty slow, so we do it only once, unless uses clears
    % everything
    [MY_LICENSE_KEY,Err] = mxlindo('LSloadLicenseString',MY_LICENSE_FILE);
    [iEnv,nErr]=mxlindo('LScreateEnv',MY_LICENSE_KEY);
    if nErr ~= LSERR_NO_ERROR;output = returnempty(-5); return; end;
end

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

% Do some pre-calc to be used in callbacks
nonlinearindicies = find(interfacedata.variabletype~=0);
nonlinearindicies = unionstripped(nonlinearindicies,interfacedata.evalVariables);
linearindicies    = find(interfacedata.variabletype==0);
linearindicies    = setdiff1D(linearindicies,nonlinearindicies);
interfacedata.nonlinearindicies = nonlinearindicies;
interfacedata.linearindicies    = linearindicies;

% Init model size
m  = K.l + K.f;
n  = length(c);

% Specifying variable types...
vtype = repmat('C',1,length(c(linearindicies)));
vtype(interfacedata.integer_variables) = 'I';

oshift = interfacedata.f;

if m == 0
    interfacedata.F_struc = [1e6 -ones(1,length(c))];
    K.l = 1;
    F_struc = [1e6 -ones(1,length(c))];
    m = 1;
    csense = [repmat('E',1,K.f) repmat('L',1,K.l)];
end

[Nbegcol,Nlencol,Nrowndx] = lindosparse(ones(m,length(linearindicies)));
oJacobian = ones(length(linearindicies),1);
Nobjndx = find(oJacobian) - 1;
Nobjcnt = length(Nobjndx);
if  isempty(Nobjndx)
    Nobjndx = [];
end

%[Nbegcol,Nlencol,Nrowndx,Nobjcnt,Nobjndx,Apatt] = jacSparsityGeometric(interfacedata);
% A = -F_struc(:,1+linearindicies);
% b = full(F_struc(:,1));
csense = [repmat('E',1,K.f) repmat('L',1,K.l)];
% A = A.*(~Apatt);
% b(any(Apatt,2)) = 0;

[iModel,nErr]=mxlindo('LScreateModel',iEnv);
if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;

% SETUP data for callbacks
mt      = interfacedata.monomtable;
linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
nonlinear_variables = setdiff((1:size(mt,1))',linear_variables);
sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));
extended_variables = interfacedata.extended_variables;

[prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);          
prob.interfacedata = interfacedata;
if problem   
    output = createoutput([],[],[],problem,'LINDO',[],[],0);
    return
end
lindo_fungp([],[],[],[],[],[],prob);

[nErr] = mxlindo('LSsetFuncalc', iModel, 'lindo_fungp',prob);
if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;

[nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_PRINTLEVEL, options.verbose+1);
if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;

% Set NLP solver
[nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_SOLVER, eval(options.lindo.LS_IPARAM_NLP_SOLVER));
[nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_MAXLOCALSEARCH,options.lindo.LS_IPARAM_NLP_MAXLOCALSEARCH);

% Load the LP portion of  model
if ~isempty(lb)
    lb = lb(linear_variables);
    ind = find(lb<1e-2);
    lb(ind) = exp(log(1e-2)+(lb(ind)-1e-2)/1e-2);
    lb = log(lb+sqrt(eps));
end
if ~isempty(ub)
    ub = ub(linear_variables);
    ind = find(ub<1e-2);
    ub(ind) = exp(log(1e-2)+(ub(ind)-1e-2)/1e-2);
    ub = log(ub+sqrt(eps));
end
b = zeros(m,1);
A = spalloc(m,length(linearindicies),0);
[nErr] = mxlindo('LSXloadLPData', iModel, 1, 0, zeros(length(linearindicies),1),b, csense,A, lb, ub);
if nErr ~= LSERR_NO_ERROR;output = createoutput(11); return; end;

nErr = mxlindo('LSloadVarType',iModel,vtype);
if nErr ~= LSERR_NO_ERROR;output = createoutput(11); return; end;

% Load the NLP portion of the model
[nErr] = mxlindo('LSloadNLPData', iModel, Nbegcol, Nlencol,[], Nrowndx, Nobjcnt,Nobjndx,[]);
if nErr ~= LSERR_NO_ERROR;output = createoutput(11); return; end;

% Optimize model
solvertime = clock;

if isempty(interfacedata.integer_variables)
    solver = 2;
else
    solver = 1;
end
solvertime = clock;
switch solver
    case 1
        [solstat,nErr] = mxlindo('LSsolveMIP', iModel);
        if ~ismember(solstat,[2009 LS_STATUS_INFEASIBLE])
            [x,nErr] = mxlindo('LSgetMIPPrimalSolution',iModel);
        else
            x = zeros(length(linearindicies),1);
        end
    case 2
        [solstat,nErr] = mxlindo('LSoptimize', iModel, eval(options.lindo.LS_METHOD));
        if ~ismember(solstat,[2009 LS_STATUS_INFEASIBLE])
            [x,nErr] = mxlindo('LSgetPrimalSolution',iModel);
        else
            x = zeros(length(linearindicies),1);
        end
    case 3
        [solStatus,nErr] = mxlindo('LSsolveGOP', iModel);
        [x,nErr] = mxlindo('LSgetPrimalSolution',iModel);
    otherwise
end
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

w = zeros(length(c),1);w(linearindicies) = exp(x);
y = [];

%[nErr]=mxlindo('LSdeleteEnv',iEnv);
[nErr]=mxlindo('LSdeleteModel',iModel);

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
    solverinput.solstat = solstat;
    solverinput.nErr = nErr;
    solverinput.x = x;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
    solveroutput.lambda=lambda;
else
    solveroutput = [];
end

% Standard interface
output = createoutput(w,y,[],problem,'LINDO',solverinput,solveroutput,solvertime);
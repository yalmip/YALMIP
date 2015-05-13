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

% % Move nonlinear bounds to constraints
if ~isempty(lb)
    finite = find(~isinf(lb(nonlinearindicies)));
    if ~isempty(finite)
        temp = F_struc(1:K.f,:);
        F_struc(1:K.f,:) = [];
        n = length(c);

        for i = 1:length(finite)
            j = nonlinearindicies(i);
            F_struc = [-lb(j) sparse(1,j,1,1,n);F_struc];
        end
        K.l = K.l + length(finite);
        F_struc = [temp;F_struc];
        interfacedata.K = K;
        interfacedata.F_struc = F_struc;
    end
end
if ~isempty(ub)
    finite = find(~isinf(ub(nonlinearindicies)));
    if ~isempty(finite)
        temp = F_struc(1:K.f,:);
        F_struc(1:K.f,:) = [];
        n = length(c);

        for i = 1:length(finite)
            j = nonlinearindicies(i);
            F_struc = [ub(j) -sparse(1,j,1,1,n);F_struc];
        end
        K.l = K.l + length(finite);
        F_struc = [temp;F_struc];
        interfacedata.K = K;
        interfacedata.F_struc = F_struc;
    end
end
% 
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

[Nbegcol,Nlencol,Nrowndx,Nobjcnt,Nobjndx,Apatt] = jacSparsity(interfacedata);
A = -F_struc(:,1+linearindicies);
b = full(F_struc(:,1));
csense = [repmat('E',1,K.f) repmat('L',1,K.l)];
A = A.*(~Apatt);
b(any(Apatt,2)) = 0;


%[MY_LICENSE_KEY,nErr] = mxlindo('LSloadLicenseString',MY_LICENSE_FILE);

%[iEnv,nErr]=mxlindo('LScreateEnv',MY_LICENSE_KEY);
%if nErr ~= LSERR_NO_ERROR;output = returnempty(-5); return; end;
[iModel,nErr]=mxlindo('LScreateModel',iEnv);
if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
constant_data = setup_fmincon_params(interfacedata);
constant_data.F_struc = F_struc;
lindo_fun([],[],[],[],[],[],constant_data);
[nErr] = mxlindo('LSsetFuncalc', iModel, 'lindo_fun',constant_data);
if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;
[nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_PRINTLEVEL, options.verbose+1);
if nErr ~= LSERR_NO_ERROR;output = returnempty(11); return; end;

% Set NLP solver
[nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_SOLVER, eval(options.lindo.LS_IPARAM_NLP_SOLVER));
[nErr] = mxlindo('LSsetModelIntParameter', iModel, LS_IPARAM_NLP_MAXLOCALSEARCH,options.lindo.LS_IPARAM_NLP_MAXLOCALSEARCH);

% Load the LP portion of  model
[nErr] = mxlindo('LSXloadLPData', iModel, 1, 0, c(linearindicies), b, csense,sparse(A), lb(linearindicies), ub(linearindicies));
if nErr ~= LSERR_NO_ERROR;output = createoutput(11); return; end;

nErr = mxlindo('LSloadVarType',iModel,vtype);
if nErr ~= LSERR_NO_ERROR;output = createoutput(11); return; end;

% Load the NLP portion of the model
[nErr] = mxlindo('LSloadNLPData', iModel, Nbegcol, Nlencol,[], Nrowndx, Nobjcnt,Nobjndx,[]);
if nErr ~= LSERR_NO_ERROR;output = createoutput(11); return; end;

% Optimize model
solvertime = tic;

if isempty(interfacedata.integer_variables)
    solver = 2;
else
    solver = 1;
end
solvertime = tic;
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
solvertime = toc(solvertime);

w = zeros(length(c),1);w(linearindicies) =x;
y = [];

%[nErr]=mxlindo('LSdeleteEnv',iEnv);
[nErr]=mxlindo('LSdeleteModel',iModel);

switch solstat
    case {LS_STATUS_OPTIMAL,LS_STATUS_BASIC_OPTIMAL,7,8}
        problem = 0;
    case {LS_STATUS_INFEASIBLE,LS_STATUS_LOCAL_INFEASIBLE}
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
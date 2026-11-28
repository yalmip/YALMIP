function output = callglpkoctave(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
interfacedata.gettime = 0;
n = length(c);

% Remember where the artificial bounds +-1e6 below are added: a solution at one of them
% means that the problem is unbounded
artificialLB = isempty(ub);
artificialUB = isempty(ub);
if ~isempty(ub)
    LB = lb;
    UB = ub;
    LB(binary_variables)  = round(LB(binary_variables));
    UB(binary_variables)  = round(UB(binary_variables));
    LB(integer_variables) = round(LB(integer_variables));
    UB(integer_variables) = round(UB(integer_variables));

    if all(isinf(LB))
        LB=repmat(-1e6,n,1);    % just for sure
        artificialLB = true;
    end
    if all(isinf(UB))
        UB=repmat(1e6,n,1);    % just for sure
        artificialUB = true;
    end
else
    LB=repmat(-1e6,n,1);    % just for sure
    UB=repmat(1e6,n,1);
end

% GLPK does not like lb==ub
equality_in_bound = find((abs(LB-UB)<1e-12) & ~isinf(LB));
m = length(equality_in_bound);
if ~isempty(equality_in_bound)
    F_struc = [-LB(equality_in_bound) sparse(1:m,equality_in_bound,ones(m,1),m,n);F_struc];
    UB(equality_in_bound) = UB(equality_in_bound) + 1;
    LB(equality_in_bound) = LB(equality_in_bound) - 1;
    K.f = K.f + m;
end

if options.showprogress;showprogress('Calling GLPK',options.showprogress);end

% GLPK notation for sure...
SENSE = 1;     % Minimize
C = full(c);   % Must be full
B = full(F_struc(:,1));         % Must be full
A =-F_struc(:,2:end);
if length(B)==0;
    A = C';
    B = 1e6;
end
% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
CTYPE = [char(ones(K.f,1)*83); char(ones(K.l,1)*85)];
VARTYPE = char(ones(n,1)*67);
VARTYPE(integer_variables) = 'I'; 
VARTYPE(binary_variables)  = 'I';  % Should not happen except from bmibnb

if options.savedebug
    save glpkmexdebug
end
options.glpk.msglev = options.verbose;
if options.glpk.msglev==1
    options.glpk.msglev = 2;
end

% Call mex-interface
solvertime = tic;
[x,FMIN,STATUS,LAMBDA_EXTRA] = glpk(C,A,B,LB,UB,CTYPE,VARTYPE,SENSE);
solvertime = toc(solvertime);

problem = 0;

if options.saveduals
try
    if isstruct(LAMBDA_EXTRA)    
        LAMBDA = LAMBDA_EXTRA.lambda;
        EXTRA = LAMBDA_EXTRA;
    else
        LAMBDA = LAMBDA_EXTRA;
        EXTRA = 'Not saved in old GLPKMEX version, update...';
    end
    D_struc = -LAMBDA(m+1:end,:);
    catch
    D_struc = [];
    end
else
    D_struc = [];
end

% STATUS is GLPK's error code; on success, LAMBDA_EXTRA.status holds the solution status
% (GLP_OPT 5, GLP_FEAS 2, GLP_NOFEAS 4, GLP_UNBND 6)
solstatus = 5;
if isstruct(LAMBDA_EXTRA) && isfield(LAMBDA_EXTRA,'status')
    solstatus = LAMBDA_EXTRA.status;
end
% An artificial bound that is binding means that the problem is unbounded. For an LP,
% the bound is binding when its reduced cost is nonzero (positive at a lower bound,
% negative at an upper bound); otherwise fall back to the sign of the cost.
atLB = artificialLB & abs(x-LB) <= 1e-6;
atUB = artificialUB & abs(x-UB) <= 1e-6;
if isstruct(LAMBDA_EXTRA) && isfield(LAMBDA_EXTRA,'redcosts') && numel(LAMBDA_EXTRA.redcosts) == n ...
        && isempty(integer_variables) && isempty(binary_variables)
    d = LAMBDA_EXTRA.redcosts(:);
    tol = 1e-9*max(1,norm(C,inf));
    atartificial = any(atLB & d > tol) || any(atUB & d < -tol);
else
    atartificial = any(atLB & C > 0) || any(atUB & C < 0);
end

% Check, currently not exhaustive...
if STATUS == 0 && (solstatus == 6 || atartificial)
    problem = 2;        % unbounded
elseif STATUS == 0 && solstatus == 4
    problem = 1;        % infeasible
elseif STATUS == 0
    problem = 0;
elseif STATUS == 10
    problem = 1;        % presolver: no primal feasible solution
elseif STATUS == 11
    problem = 12;       % presolver: no dual feasible solution, infeasible or unbounded
else
    problem = -1;
end

% Save all data sent to solver?
if options.savesolverinput
	solverinput.A = A;
	solverinput.C = C;
	solverinput.B = B;
	solverinput.CTYPE = CTYPE;
	solverinput.LB = LB;
	solverinput.UB = UB;
    solverinput.param = options.glpk;
    solverinput.lpsolver = options.glpk.lpsolver;
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput

    % Stupid redo, but this saves some time in bnb calls
    if isstruct(LAMBDA_EXTRA)
        LAMBDA = LAMBDA_EXTRA.lambda;
        EXTRA = LAMBDA_EXTRA;
    else
        LAMBDA = LAMBDA_EXTRA;
        EXTRA = 'Not saved in old GLPKMEX version, update...';
    end
    
	solveroutput.x = x;
    solveroutput.FMIN = FMIN;
    solveroutput.STATUS = STATUS;
    solveroutput.LAMBDA=LAMBDA;
    solveroutput.EXTRA = EXTRA; 
else
	solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,yalmiperror(problem,interfacedata.solver.tag),solverinput,solveroutput,solvertime);
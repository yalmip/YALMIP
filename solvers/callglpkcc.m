function output = callglpk(interfacedata)

% Author Johan Löfberg 
% $Id: callglpkcc.m,v 1.2 2008-03-28 12:04:51 joloef Exp $

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

if ~isempty(ub)
    LB = lb;
    UB = ub;
    LB(binary_variables)  = round(LB(binary_variables));
    UB(binary_variables)  = round(UB(binary_variables));
    LB(integer_variables) = round(LB(integer_variables));
    UB(integer_variables) = round(UB(integer_variables));

    if all(isinf(LB))
        LB=repmat(-1e6,n,1);    % just for sure
    end
    if all(isinf(UB))
        UB=repmat(1e6,n,1);    % just for sure
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
    K.l = 1;
end
% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
CTYPE = [char(ones(K.f,1)*83); char(ones(K.l,1)*85)];
VARTYPE = char(ones(n,1)*67);
VARTYPE(integer_variables) = 'I'; 
VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb

if options.savedebug
    save glpkmexdebug
end
options.glpk.msglev = options.verbose;
if options.glpk.msglev==1
    options.glpk.msglev = 2;
end

% Call mex-interface
solvertime = clock; 
%[x,FMIN,STATUS,LAMBDA_EXTRA] = glpkmex(SENSE,C,A,B,CTYPE,LB,UB,VARTYPE,options.glpk,options.glpk.lpsolver,options.glpk.save);
[x,FMIN,STATUS,LAMBDA_EXTRA] = glpkcc(C, A, B, LB, UB, CTYPE,VARTYPE,SENSE,options.glpk);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end
problem = 0;

if options.saveduals
    if isstruct(LAMBDA_EXTRA)
        LAMBDA = LAMBDA_EXTRA.lambda;
        EXTRA = LAMBDA_EXTRA;
    else
        LAMBDA = LAMBDA_EXTRA;
        EXTRA = 'Not saved in old GLPKMEX version, update...';
    end
    D_struc = -LAMBDA(m+1:end,:);
else
    D_struc = [];
end

% Hack
if isempty([ub;lb]) & (any(LB==x & C>0) | any(UB==x & C<0))
    STATUS = 214;
end

% Check, currently not exhaustive...
switch STATUS
    case {5}
        problem = 0;
    case {3,4,110,213}
        problem = 1;
    case {6,214}
        problem = 2;
    case 207
        problem = 3;
    case {210,211,212}
        problem = 4;
    case {1,170}
        problem = 11;
    otherwise
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
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;
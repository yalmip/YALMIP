function output = callscipmex(interfacedata)

% Author Johan Löfberg 

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

if options.showprogress;showprogress('Calling SCIP',options.showprogress);end

f = full(c);   % Must be full
ru = full(F_struc(:,1));         % Must be full
A =-F_struc(:,2:end);
rl = -inf(length(ru),1);
rl(1:K.f) = ru(1:K.f);

% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
VARTYPE = char(ones(n,1)*67);
VARTYPE(integer_variables) = 'I'; 
VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb

if options.savedebug
    save scipmexdebug
end

% Call mex-interface
solvertime = clock; 
%[x,FMIN,STATUS,LAMBDA_EXTRA] = glpkmex(SENSE,C,A,B,CTYPE,LB,UB,VARTYPE,options.glpk,options.glpk.lpsolver,options.glpk.save);
[x,FMIN,STATUS,INFO] = scip(f, A, rl,ru, lb, ub,VARTYPE);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end
problem = 0;

D_struc = [];

% Check, currently not exhaustive...
switch STATUS
    case {2,3,4,5}
        problem = 3;
    case {10}
        problem = 0;   
    case 11
        problem = 1;
    case 12
        problem = 2;
    case 13
        problem = 15;
    case {1,6,7,8,9}
        problem = 11;
    otherwise
        problem = -1;
end

% Save all data sent to solver?
if options.savesolverinput
	solverinput.A = A;
	solverinput.f = f;
	solverinput.rl = rl;
	solverinput.ru = ru;
	solverinput.lb = lb;
	solverinput.ub = ub;
    solverinput.VARTYPE = VARTYPE;
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
	solveroutput.x = x;
    solveroutput.FMIN = FMIN;
    solveroutput.STATUS = STATUS; 
    solveroutput.INFO = INFO;   
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
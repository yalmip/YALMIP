function output = callglpk(interfacedata)

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
    lb(binary_variables)  = round(lb(binary_variables));
    ub(binary_variables)  = round(ub(binary_variables));
    lb(integer_variables) = round(lb(integer_variables));
    ub(integer_variables) = round(ub(integer_variables));

    if all(isinf(lb))
        lb=repmat(-1e6,n,1);    % just for sure
    end
    if all(isinf(ub))
        ub=repmat(1e6,n,1);    % just for sure
    end
else
  lb=repmat(-1e6,n,1);   
  ub=repmat(1e6,n,1);    
end

if options.showprogress;showprogress('Calling OSL',options.showprogress);end

B = full(F_struc(K.f+1:end,1));         
A =-F_struc(K.f + 1:end,2:end);
Beq = full(F_struc(1:K.f:end,1));         
Aeq =-F_struc(1:K.f:end,2:end);
if length(B)==0;
    A = C';
    B = 1e6;
end

if options.savedebug
    save osldebug
end

STATUS = 0;
solvertime = tic;
[x,fval] = oslprog(c,A,B,Aeq,Beq,lb,ub);
solvertime = toc(solvertime);
problem = 0;

% No duals returned
D_struc = [];

% Check, currently not exhaustive...
switch STATUS
    case 0
        problem = 0;
    otherwise
        problem = -1;
end

% Save all data sent to solver?
if options.savesolverinput
	solverinput.A = A;
	solverinput.C = C;
	solverinput.B = B;
	solverinput.Aeq = Aeq;
    solverinput.beq = beq;
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

if ~options.savesolverinput
    solverinput = [];
else
    solverinput = model;
end
if ~options.savesolveroutput
    solveroutput = [];
else
    solveroutput = solveroutput;
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);
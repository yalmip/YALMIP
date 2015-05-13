function output = callcplexmex(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
LB      = interfacedata.lb;
UB      = interfacedata.ub;

showprogress('Calling CPLEXMEX',options.showprogress);

if ~isempty(LB)
    LB(isinf(LB)) = -1e12;
    UB(isinf(UB)) = 1e12;
end

SENSE = 1;   
C = c(:);
if ~isempty(F_struc)
    A =-F_struc(:,2:end);
    B = full(F_struc(:,1));
else
    A = zeros(1,length(c));A(1)=1;
    B = 1e6;
end

H = 2*full(Q);
if nnz(H)==0
    H = [];
end

CTYPE = repmat('L',K.l+K.f,1);  % Standard variables
CTYPE(1:K.f) = 'E';             % Equality constrained variables
VARTYPE = repmat('C',size(A,2),1);
VARTYPE(binary_variables)  = 'B';
VARTYPE(integer_variables) = 'I';

options.cplexmex.msglev = options.verbose;
if options.cplexmex.msglev>1
    options.cplexmex.msglev = 1;
end

% Call mex-interface
if options.savedebug
save cplexmexdebug
end
solvertime = tic;
[x,OPT,STATUS,EXTRA]= cplexmex(SENSE,H,C,A,B,CTYPE,LB,UB,VARTYPE,x0,options.cplexmex);
solvertime = toc(solvertime);
problem = 0;

D_struc = -EXTRA.lambda;

% Check, currently not exhaustive...
switch STATUS
    case {1,101,102}
        problem = 0;
    case {3,103}
        problem = 1;
    case {2,118}
        problem = 2;
    case {4,119}
        problem = 12;
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'CPLEXMEX');	

% Save all data sent to solver?
if options.savesolverinput
    solverinput.SENSE = SENSE;    
    solverinput.H = H;
	solverinput.A = A;
	solverinput.C = C;
	solverinput.B = B;
	solverinput.CTYPE = CTYPE;
	solverinput.LB = LB;
	solverinput.UB = UB;
	solverinput.VARTYPE = VARTYPE;
    solverinput.x0 = x0;    
    solverinput.param = options.cplexmex;
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
	solveroutput.X = x;
    solveroutput.OPT = OPT;
    solveroutput.STATUS = STATUS;
    solveroutput.EXTRA = EXTRA; 
else
	solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);
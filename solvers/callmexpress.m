function output = callmexpress(interfacedata)

% Author Johan Löfberg 
% $Id: callmexpress.m,v 1.7 2005-05-07 13:53:20 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
UB      = interfacedata.ub;
LB      = interfacedata.lb;

showprogress('Calling XPRESS',options.showprogress);

SENSE = 1;     
C = full(c);   
if isempty(F_struc)
    A = zeros(1,length(c));A(1)=1;
    B = 1e6;
else
    A =-(F_struc(:,2:end)); 
    B = full(F_struc(:,1));            
end

% XPRESS complains about large bounds
if isempty(LB)
  %  LB = repmat(-2147483647,length(c),1);
else
    LB(LB<-2147483647) = -2147483647;
end
if isempty(UB)
  %  UB = repmat(2147483647,length(c),1);
else
    UB(UB>2147483647) = 2147483647;
end

CTYPE = repmat('L',K.l+K.f,1);  % Inequalities
CTYPE(1:K.f) = 'E';             % Equality constraints
VARTYPE = repmat('C',size(A,2),1);
VARTYPE(integer_variables)='I'; % Integer variables
VARTYPE(binary_variables) ='B'; % Binary variables

H = full(2*Q);
if nnz(Q)==0
    H = [];
end

if options.verbose==0
    options.xpress.msglev = 0;
else
    options.xpress.msglev = 1;
end

if options.savedebug
    save mexpressdebug SENSE H C A B CTYPE LB UB VARTYPE options
end

% Call mex-interface
solvertime = clock; 
[x,FMIN,STATUS,EXTRA] = mexpress(SENSE,H,C,A,B,CTYPE,LB,UB,VARTYPE,options.xpress,0);
solvertime = etime(clock,solvertime);
problem = 0;
if isstruct(EXTRA)
    D_struc = -EXTRA.lambda;    
else
    D_struc = [];
end

% Check (error code depends on problem type!)
if isempty(union(binary_variables,integer_variables))
    switch STATUS
        case 1
            problem = 0;
        case 2
            problem = 1;
        case 5
            problem = 2;
        case 4
            problem = 3;
        case {0,3,6}
            problem = 11;
        otherwise
            problem = -1;
    end
else
    switch STATUS
        case 6
            problem = 0;
        case 5
            problem = 1;
        case 4
            problem = 2;
        case {4,3}
            problem = 3;
        case {0,1,2}
            problem = 11;
        otherwise
            problem = -1;
    end
    
end
infostr = yalmiperror(problem,'MEXPRESS');	

% Save all data sent to solver?
if options.savesolverinput
	solverinput.H = H;
    solverinput.A = A;
	solverinput.C = C;
	solverinput.B = B;
	solverinput.CTYPE = CTYPE;
	solverinput.LB = LB;
	solverinput.UB = UB;
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
	solveroutput.x = x;
    solveroutput.FMIN = FMIN;
    solveroutput.STATUS = STATUS;
    solveroutput.EXTRA=EXTRA;
else
	solveroutput = [];
end


% Standard interface 
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;
function output = callscipmex(interfacedata)

% Fix for the case when YALMIP selects scip as a solver for a nonlinearly
% parameterized problem, and te final problem is actually still nonlonear.
% If so, we should really have selected scipnl
if any(interfacedata.variabletype) || ~isempty(interfacedata.evalMap)
    output = callscipnl(interfacedata);
    return
end

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
H       = 2*interfacedata.Q;
interfacedata.gettime = 0;
n = length(c);

if options.showprogress;showprogress('Calling SCIP',options.showprogress);end

f = full(c);   % Must be full
if K.f+K.l == 0
    A = [];
    ru = [];
    rl = [];
else
    ru = full(F_struc(1:K.f+K.l,1));         % Must be full
    A =-F_struc(1:K.f+K.l,2:end);
    rl = -inf(length(ru),1);
    rl(1:K.f) = ru(1:K.f);
end

% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
VARTYPE = char(ones(n,1)*67);
VARTYPE(integer_variables) = 'I'; 
VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb

if any(K.q)
    top = K.f+K.l+1; 
    qc.l = [];
    qc.r = [];
    for i = 1:length(K.q)
        ni = K.q(i);
        % |Qx+d|<a+b'x
        ai = F_struc(top,1);
        bi = F_struc(top,2:end);
        Qi = F_struc(top+1:top+K.q(i)-1,2:end);
        di = F_struc(top+1:top+K.q(i)-1,1);
        qc.Q{i} = Qi'*Qi-bi'*bi;
        qc.l = [qc.l full(2*Qi'*di-2*ai*bi')];
        qc.qru = [qc.r;full(ai*ai-di'*di)];
        qc.qrl = [qc.r;-50];
        % Linear constraint a+b'*x >= 0 should be added
        top = top  + K.q(i);
    end
    if i == 1
        qc.Q = qc.Q{1};
    end
else
    qc = [];
end

sos.type = K.sos.type;
if isempty(sos.type)
    sos.type = '';
end
sos.index = K.sos.variables;
sos.weight = K.sos.weight;

ops = options.scip;
switch options.verbose
    case 0
        ops.display = 0;
    case 3
        ops.display = 3;
    otherwise
        ops.display = 4;
end

if options.savedebug
    save scipdebug H f A rl ru lb ub VARTYPE sos ops
end

% Call mex-interface
solvertime = tic;
try
    [x,FMIN,STATUS,INFO] = scip(H, f, A, rl, ru, lb, ub, VARTYPE, sos,qc,[],ops);
catch
    % -6 for instance causes a blank crash
    x = [];
    FMIN=[];
    STATUS = -1;
    INFO = [];
end
solvertime = toc(solvertime);

D_struc = [];

% Check, currently not exhaustive...
problem = 0;
switch STATUS
    case 0
        problem = 9;
    case 1
        problem = 16;
    case {2,3,4,5}
        problem = 3;
    case {7,8,9,10}
        problem = 4;
    case 11
        problem = 0;
    case 12
        problem = 1;
    case 13
        problem = 2;
    case 14
        problem = 15;
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
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
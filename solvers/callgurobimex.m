function output = callgurobi(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
Q       = interfacedata.Q;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
x0      = interfacedata.x0;
interfacedata.gettime = 0;
n = length(c);

if ~isempty(ub)
    LB = lb;
    UB = ub;
    LB(binary_variables)  = round(LB(binary_variables));
    UB(binary_variables)  = round(UB(binary_variables));
    LB(integer_variables) = round(LB(integer_variables));
    UB(integer_variables) = round(UB(integer_variables));
else
    LB = [];
    UB = [];
end

if options.showprogress;showprogress('Calling GUROBI',options.showprogress);end

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [LB,UB,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,LB,UB);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
    
    redundant = find(LB<=0 & UB>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
        
end

SENSE = 1;     % Minimize
C = full(c);   % Must be full
if size(F_struc,1)>0
    B = full(F_struc(:,1));         % Must be full
    A =-F_struc(:,2:end);
else
    B = [];
    A = [];
end

% Optimized code, make a lot of difference when you make this call 10000
% times in a branch and bound setting...
CTYPE = [char(ones(K.f,1)*61); char(ones(K.l,1)*60)];
VARTYPE = char(ones(length(c),1)*67);
VARTYPE(setdiff(integer_variables,semicont_variables)) = 'I';
VARTYPE(binary_variables)  = 'B';  % Should not happen except from bmibnb
VARTYPE(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
VARTYPE(intersect(semicont_variables,integer_variables)) = 'N';

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(UB(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = UB(semicont_variables(NegativeSemiVar));
        UB(semicont_variables(NegativeSemiVar)) = -LB(semicont_variables(NegativeSemiVar));
        LB(semicont_variables(NegativeSemiVar)) = -temp;
        A(:,semicont_variables(NegativeSemiVar)) = -A(:,semicont_variables(NegativeSemiVar));
        C(semicont_variables(NegativeSemiVar)) = -C(semicont_variables(NegativeSemiVar));
        if ~isempty(x0)
            x0(NegativeSemiVar) = -NegativeSemiVar;
        end
    end
end

if nnz(Q)>0
    [ii,jj,kk] = find(Q);
    ii = ii-1;
    jj = jj-1;
    comp = computer;
    % According to Wotao's testing
    if isempty(strfind(comp,'64')) | isequal(comp,'GLNXA64') | isequal(comp,'PCWIN64')
        options.gurobi.QP.qrow = int32(ii)';
        options.gurobi.QP.qcol = int32(jj)';
    elseif strfind(comp,'MACI64')
        options.gurobi.QP.qrow = int32(ii)';
        options.gurobi.QP.qcol = int32(jj)';    
    else
        options.gurobi.QP.qrow = int64(ii)';
        options.gurobi.QP.qcol = int64(jj)';
    end
    options.gurobi.QP.qval = kk';
end
        
if ~options.verbose
    options.gurobi.DisplayInterval = 0;
    options.gurobi.Display = 0;
end

if ~isempty(x0)
    options.gurobi.Start = x0(:)';
end

if isinf(options.gurobi.SolutionLimit)
    options.gurobi.SolutionLimit = double(intmax);
end
if isinf(options.gurobi.BarIterLimit)
    options.gurobi.BarIterLimit = double(intmax);
end
if isfield(options.gurobi,'InfUnbdInfo')
    % Options really for Gurobi 5.0 
    options.gurobi=rmfield(options.gurobi,'InfUnbdInfo');
end

if options.savedebug
    save gurobidebug
end

if ~isempty(K.sos.type)
    options.gurobi.SOS.weights = spalloc(length(c),length(K.sos.type),0);
    for i = 1:length(K.sos.type)
        options.gurobi.SOS.types(i)= int32(str2num(K.sos.type(i)));
        options.gurobi.SOS.weights(K.sos.variables{i},i) = 1:length(K.sos.variables{i});
    end
end

% Call mex-interface
solvertime = tic;
if isempty(binary_variables) & isempty(integer_variables) & isempty(semicont_variables)
    [x,val,flag,output,lambda] = gurobi_mex(C,SENSE,sparse(A),B,CTYPE,LB,UB,VARTYPE,options.gurobi);
else
    [x,val,flag,output] = gurobi_mex(C,SENSE,sparse(A),B,CTYPE,LB,UB,VARTYPE,options.gurobi);
    lambda = [];
end
solvertime = toc(solvertime);

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
if length(x) == length(C)
    if ~isempty(NegativeSemiVar)
        x(NegativeSemiVar) = -x(NegativeSemiVar);
    end
end

problem = 0;
if isa(lambda,'double')
    D_struc = -lambda;
elseif isa(lambda,'struct')
    D_struc = -lambda.Pi;
end

if isempty(flag)
    flag = 9;
    x = nan(length(C),1);
end

% Check, currently not exhaustive...
switch flag
    case {2}
        problem = 0;
    case {3}
        x = zeros(length(C),1);
        problem = 1;
    case 4
        x = zeros(length(C),1);
        problem = 12;
    case 5
        problem = 2;
    case {7,8,9}
        problem = 3;
    case {12,13}
        problem = 4;
    case {1,6,10,11}
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
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput    
	solveroutput.x = x;
    solveroutput.val = val;
    solveroutput.flag = flag;
    solveroutput.lambda=lambda;
    solveroutput.output = output; 
else
	solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);
function output = calllpsolvenative(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

n = length(c);
% Bounded variables converted to constraints
if ~isempty(ub)
    LB = lb;
    UB = ub;
else
    if isempty(integer_variables) & isempty(binary_variables)         
        LB = -ones(n,1)*inf;;
        UB = ones(n,1)*inf;;
    else
        %LP_SOLVE FAILS IF BOUNDS NOT EXPLICIT
        [LB,UB,used_rows] = find_lp_bounds(F_struc,K);
        LB(isinf(LB)) = -1e12;
        UB(isinf(UB)) = 1e12;
        F_struc(K.f+used_rows,:)=[];
        K.l = K.l - length(used_rows);
        LB(binary_variables) = max(LB(binary_variables),0);
        UB(binary_variables) = min(UB(binary_variables),0);
    end
end

if options.showprogress;showprogress('Calling LPSOLVE',options.showprogress);end

f = - full(c);            % Must be full
A = - F_struc(:,2:end);
b = full(F_struc(:,1));   % Must be full
e = -ones(size(A,1),1);
e(1:K.f) = 0;

xint = uniquestripped([integer_variables binary_variables]);

% Call mex-interface
solvertime = clock; 

if K.f>0
    Aeq = A(1:K.f,:);
    beq = b(1:K.f);
    Alin = A(K.f+1:end,:);
    blin = b(K.f+1:end);
    [LB,UB,Alin,blin] = remove_bounds_from_Ab(Alin,blin,LB,UB);
    A = [Aeq;Alin];
    b = [beq;blin];
    e = -ones(size(A,1),1);
    e(1:K.f) = 0;
    options.saveduals = 0;
else
    [LB,UB,A,b] = remove_bounds_from_Ab(A,b,LB,UB);
    e = -ones(size(A,1),1);
    options.saveduals = 0;
end
if ~isempty(semicont_variables)
    redundant = find(LB<=0 & UB>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
end

% LPSOLVE assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(UB(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = UB(semicont_variables(NegativeSemiVar));
        UB(semicont_variables(NegativeSemiVar)) = -LB(semicont_variables(NegativeSemiVar));
        LB(semicont_variables(NegativeSemiVar)) = -temp;
        A(:,semicont_variables(NegativeSemiVar)) = -A(:,semicont_variables(NegativeSemiVar));
        f(semicont_variables(NegativeSemiVar)) = -f(semicont_variables(NegativeSemiVar));
    end
end
if options.savedebug
    save lpsolvedebug f A b e UB LB xint 
end

lp = create_lp_solve_native_model(A,b,f,xint,LB,UB,e,options);
if ~isempty(K.sos)
    for i = 1:length(K.sos.type)
       mxlpsolve('add_SOS', lp, 'Dummy', str2num(K.sos.type(i)), 1, K.sos.variables{i}, K.sos.weight{i});
    end
end
if ~isempty(semicont_variables)
    mxlpsolve('set_semicont', lp, semicont_variables) 
end

try    
    if options.showprogress;showprogress(['Calling ' interfacedata.solver.tag],options.showprogress);end
    solvertime = tic;
    result=mxlpsolve('solve', lp);
    solvertime = toc(solvertime);
    if result == 0 | result == 1 | result == 11 | result == 12        
        [obj, x, duals] = mxlpsolve('get_solution', lp);
    else
        obj = [];
        x = zeros(length(c),1);
        duals = [];
    end
    mxlpsolve('delete_lp', lp);
catch
    obj = [];
    x = zeros(length(c),1);
    duals = [];
    result = -1;
    mxlpsolve('delete_lp', lp);
end


% LPSOLVE assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
 if length(x) == length(c)
     if ~isempty(NegativeSemiVar)
        x(NegativeSemiVar) = -x(NegativeSemiVar);
     end
 end

if options.saveduals & isempty(integer_variables)
    D_struc = duals;
else
    D_struc = [];
end

switch result
    case {0,1}
        problem = 0; % OPTIMAL
    case 2
        problem = 1; % INFEASIBLE
    case 3
        problem = 2; % UNBOUNDED
    case {7,12,13}
        problem = 3; % RUN OUT OF TIME OR SIMILIAR
    case 5
        problem = 4;
    case {-2,10,11}
        problem = 11;
    otherwise
        problem = -1;
end
        
% Save all data sent to solver?
if options.savesolverinput
	solverinput.A = A;
	solverinput.f = f;
	solverinput.b = b;	
	solverinput.LB = LB;
	solverinput.UB = UB;    
    solverinput.xint = xint;    
else
	solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
	solveroutput.x = x;
    solveroutput.obj = obj;
    solveroutput.duals = duals;  
    solveroutput.result = result;  
else
	solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
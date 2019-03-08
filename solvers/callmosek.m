function output = callmosek(model)

% Retrieve needed data
options = model.options;
F_struc = model.F_struc;
c       = model.c;
Q       = model.Q;
K       = model.K;
x0      = model.x0;
integer_variables = model.integer_variables;
binary_variables = model.binary_variables;
extended_variables = model.extended_variables;
ub      = model.ub;
lb      = model.lb;
mt      = model.monomtable;

% *********************************
% What type of variables do we have
% *********************************
model.linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
model.nonlinear_variables = setdiff((1:size(mt,1))',model.linear_variables);
model.sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));

% Some meta-solver thinks we handle binaries
if ~isempty(model.binary_variables)
    integer_variables = union(model.integer_variables, model.binary_variables);
    if isempty(lb)
        lb = repmat(-inf,1,length(model.c));
    end
    lb(model.binary_variables) = max(lb(model.binary_variables),0);    
    if isempty(ub)
        ub = repmat(inf,1,length(model.c));
    end
    ub(model.binary_variables) = min(ub(model.binary_variables),1);
    model.lb = lb;
    model.ub = ub;
    model.integer_variables = integer_variables;
end

% Some meta solvers might construct model with empty cones
if any(model.K.s) && any(model.K.s == 0)
    model.K.s(model.K.s==0)=[];
end
if any(model.K.q) && any(model.K.q == 0)
    model.K.q(model.K.q==0)=[];
end

if ~isempty(model.sigmonial_variables) | isequal(model.solver.version,'GEOMETRIC')
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_geometric(model);          
else
    [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdp(model);        
end

% YALMIP has introduced internal variables for socp/exp cones etc
if length(x) > 0 && length(x) ~= length(model.c)
    x = x(1:length(model.c));
end

infostr = yalmiperror(problem,'MOSEK');	

% Save all data sent to solver?
if options.savesolverinput
    solverinput.prob = prob;
    solverinput.param = options.mosek;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.r = r;
    solveroutput.res = res;
    
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);


function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_lpqpsocpsdp(model);

[model,output] = normalizeExponentialCone(model);
if output.problem
    problem = output.problem;
    x = [];
    D_struc = [];
    r = [];
    res = [];
    solvertime = 0;
    prob = [];
else
    if nnz(model.Q)==0 && isempty(model.integer_variables) && isempty(model.x0)
        % Standard cone problem which we can model by sending our standard dual
        % and then recover solution via Moseks dual
        [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_dual(model);
    else
        % Integer conic program
        % Quadratic objective
        % Exponential cones
        [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_primal(model);
    end
end

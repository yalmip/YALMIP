function output = callgurobi(interfacedata)

options = interfacedata.options;
nOriginal = length(interfacedata.c);
model = yalmip2gurobi(interfacedata);

if interfacedata.options.savedebug
    save gurobidebug model
end

if options.showprogress;showprogress('Calling GUROBI',options.showprogress);end
solvertime = tic;
result = gurobi(model,model.params);
solvertime = toc(solvertime);

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
if isfield(result,'x')
    x = result.x;
    if length(x) == length(model.obj)
        if ~isempty(model.NegativeSemiVar)
            x(model.NegativeSemiVar) = -x(model.NegativeSemiVar);
        end
    end
    x = x(1:nOriginal);
else
    x = zeros(nOriginal,1);
end

problem = 0;
qcDual = [];
if isfield(result,'pi')
    % Gurobi has reversed sign-convention
    D_struc = -result.pi;
    if sum(interfacedata.K.q) > 0
        % YALMIP has appended equalities to model normalized SOCP variables
        % These duals should be moved to the assumed SOCP dual position
        % order: socp, equality, linear element
        % order: equality, linear element, socp
        m = sum(interfacedata.K.q);
        socpDuals = D_struc(1:m);
        D_struc(1:m)=[];
        D_struc = [D_struc(1:interfacedata.K.f+interfacedata.K.l);socpDuals];        
    end
else
    D_struc = [];
end
if isfield(result,'qcpi')
    qcDual = result.qcpi;
end

switch result.status
    case 'OPTIMAL'
        problem =0;
    case 'INFEASIBLE'
        problem = 1;
    case 'UNBOUNDED'
        problem = 2;
    case {'ITERATION_LIMIT','TIME_LIMIT','NODE_LIMIT'}
        problem = 3;
    case {'NUMERIC','SUBOPTIMAL'}
        problem = 4;
    case 'INF_OR_UNBD'
        problem = 12;
    otherwise
        problem = -1;
end

% Save all data sent to solver?
if interfacedata.options.savesolverinput
	solverinput.model = model;
else
	solverinput = [];
end

% Save all data from the solver?
if interfacedata.options.savesolveroutput    
	solveroutput.result = result; 
else
	solveroutput = [];
end

infostr = yalmiperror(problem,interfacedata.solver.tag);

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);
output.qcDual      = qcDual;










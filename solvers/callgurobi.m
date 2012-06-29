function output = callgurobi(interfacedata)

% Author Johan Löfberg 
% $Id: callgurobimex.m,v 1.2 2010-03-23 11:11:46 joloef Exp $

options = interfacedata.options;
n_original = length(interfacedata.c);
model = yalmip2gurobi(interfacedata);

if interfacedata.options.savedebug
    save gurobidebug model
end

if options.showprogress;showprogress('Calling GUROBI',options.showprogress);end
solvertime = clock; 
result = gurobi(model,model.params);
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

if isfield(result,'x')
    x = result.x;
    x = x(1:n_original);
else
    x = zeros(n_original,1);
end

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
if length(x) == length(model.obj)
     if ~isempty(model.NegativeSemiVar)
        x(model.NegativeSemiVar) = -x(model.NegativeSemiVar);
     end
end

problem = 0;
if isfield(result,'pi')
    if length(model.obj) == n_original
        D_struc = -result.pi;
    else
    D_struc = [];
    end
else
    D_struc = [];
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
        problem = 9;
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

% Standard interface 
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;











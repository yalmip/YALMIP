function output = callconeprog(interfacedata)

options = interfacedata.options;
model = yalmip2coneprog(interfacedata);

showprogress('Calling CONEPROG',options.showprogress);

if options.savedebug    
    save model
end

solvertime = tic;
[x,fmin,flag,output,lambda] = coneprog(model.f, model.socConstraints,model.A, model.b, model.Aeq, model.beq, model.lb, model.ub,model.ops);
solvertime = toc(solvertime);
problem = 0;

% No duals for now
D_struc = [];
if isempty(x)
    x = zeros(length(model.f),1);
end

% Check, currently not exhaustive...
if flag==0
    problem = 3;
elseif flag == -2 | flag==-5
    problem = 1;
elseif flag == -3
    problem = 2;
else
    if flag>0
        problem = 0;
    else 
        if any((model.A*x-model.b)>sqrt(eps)) | any( abs(model.Aeq*x-model.beq)>sqrt(eps))
            problem = 1; % Likely to be infeasible
        else
            if model.f'*x<-1e10 % Likely unbounded
                problem = 2;
            else          % Probably convergence issues
                problem = 5;
            end
        end
    end
end
    
% Save all data sent to solver?
if options.savesolverinput
    solverinput.model = model;   
    solverinput.options = options.coneprog;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fmin = fmin;
    solveroutput.flag = flag;
    solveroutput.output=output;
else
    solveroutput = [];
end

% Standard interface 
output = createOutputStructure(x(:),D_struc,[],problem,interfacedata.solver.tag,solverinput,solveroutput,solvertime);
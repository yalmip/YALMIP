function output = call_cplexibm_miqp(interfacedata)

%Turn on support for nonconvex QP if required and user hasn't touched this
if interfacedata.ProblemClass.objective.quadratic.nonconvex
    if isfield(interfacedata.options.cplex,'solutiontarget')
        % User is using old version of cplex
        if ~interfacedata.options.cplex.solutiontarget          
            % User has not set it to 1 or 2 manually
            interfacedata.options.cplex.solutiontarget = 3;
        end
    elseif isfield(interfacedata.options.cplex,'optimalitytarget')
        if ~interfacedata.options.cplex.optimalitytarget          
            % User has not set it to 1 or 2 manually
            interfacedata.options.cplex.optimalitytarget = 3;
        end
    else
        interfacedata.options.cplex.optimalitytarget = 3;
    end
end

options = interfacedata.options;
[model,nonlinearremain] = yalmip2cplex(interfacedata);

if nonlinearremain
    error('Nonlinear monomials remain when calling CPLEX. If you are using OPTIMIZER, ensure your model really is solvable by CPLEX for fixed parameters. If you still think so, please report this and ask for a feature improvement.');
end
if interfacedata.options.savedebug
    save cplexdebug model
end

solvertime = tic;
[x,fval,exitflag,output,lambda] = localSolverCall(model);
solvertime = toc(solvertime);
if output.cplexstatus == 4 | output.cplexstatus == 119
    % CPLEX reports infeasible OR unbounded
    % Remove objective and resolve
    model.H = model.H*0;
    model.f = model.f*0;
    solvertime = tic;
    [x,fval,exitflag,output,lambda] = localSolverCall(model);
    solvertime = toc(solvertime);
    switch output.cplexstatus
        case {1,101,102} % It was ok, hence it must have been unbounded
            output.cplexstatus = 2;
        case {3,103,106} % Infeasible, so original was infeasible
            output.cplexstatus = 3; % Infeasible
        otherwise
            output.cplexstatus = 4; % I give up
    end
end

if ~isempty(lambda)
    D_struc = [lambda.eqlin;lambda.ineqlin];
else
    D_struc = [];    
end

if size(x,1) == length(model.f)
    if ~isempty(model.NegativeSemiVar)
        x(model.NegativeSemiVar,:) = -x(model.NegativeSemiVar,:);
    end
end

showprogress('Calling CPLEX-IBM',options.showprogress);
% Check, currently not exhaustive...
switch output.cplexstatus
    case {1,101,102}
        problem = 0;
    case {3,103,106}
        problem = 1; % Infeasible
    case {2,20,21,118,133}
        problem = 2; % Unbounded
    case {4,119}
        problem = 15;    
    case {10,11,104,105,107,108,111,112}
        problem = 3; % Iteration/time
    case 119
        problem = 15;
    case {5,6,109,110}
        problem = 4; % Numerics
    case 24
        problem = 11;
    otherwise
        problem = -1;
end

infostr = yalmiperror(problem,'CPLEX-IBM');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.model = model;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.lambda = lambda;
    solveroutput.output = output;
else
    solveroutput = [];
end

if isempty(x)
    x = zeros(length(model.f),1);
end

% Standard interface 
output = createOutputStructure(x,D_struc,[],problem,infostr,solverinput,solveroutput,solvertime);

function [x,fval,exitflag,output,lambda] = localSolverCall(model)

H = model.H;
f = model.f;
Aineq = model.Aineq;
bineq = model.bineq;
Aeq = model.Aeq;
beq = model.beq;
lb = model.lb;
ub = model.ub;
x0 = model.x0;
options.cplex = model.options;
options.verbose = model.verbose;
integer_variables = model.integer_variables;
binary_variables = model.binary_variables;
semicont_variables = model.semicont_variables;
K = model.K;
ctype = model.ctype;

fixedAineqBug = 0;
if isempty(Aineq) & isempty(Aeq)
    Aineq = zeros(1,length(f)) ;
    bineq = 1;
    fixedAineqBug = 1;
end
if isempty(integer_variables) & isempty(binary_variables) & isempty(semicont_variables) & isempty(K.sos.type)
    if options.verbose 
        if isempty(H)           
            [x,fval,exitflag,output,lambda] = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);
        else
            [x,fval,exitflag,output,lambda] = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);
        end
    else
        if isempty(H)
            evalc('[x,fval,exitflag,output,lambda] = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);');            
        else
            evalc('[x,fval,exitflag,output,lambda] = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);');
        end
    end
    if ~isempty(lambda) & fixedAineqBug
        lambda.ineqlin = [];
    end
else
    if options.verbose 
        if isempty(H)
            [x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype,x0,options.cplex);
        else
            [x,fval,exitflag,output] = cplexmiqp(H,f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype,x0,options.cplex);
        end
    else
        if isempty(H)
            evalc('[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype,x0,options.cplex);');            
        else
            evalc('[x,fval,exitflag,output] = cplexmiqp(H,f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype,x0,options.cplex);');
        end
    end
    lambda = [];
end

function [x,fval,exitflag,output] = cplexmilpInternal(f,Aineq,bineq,Aeq,beq,sostype,sosvariables,sosweight,lb,ub,ctype,x0,options)

cplex = Cplex;
cplex.Model.sense = 'minimize';
cplex.Model.obj = f;
cplex.Model.lb = lb;
cplex.Model.ub = ub;
cplex.Model.ctype = ctype;
cplex.Model.A = [Aeq;Aineq];
cplex.Model.lhs = [beq;-inf(length(bineq),1)];
cplex.Model.rhs = [beq;bineq];
if options.mip.pool.intensity
    cplex.Param.mip.pool.intensity.Cur =  options.mip.pool.intensity;
    cplex.Param.mip.pool.capacity.Cur = options.mip.pool.capacity;
    cplex.Param.mip.pool.absgap.Cur = options.mip.pool.absgap;
    cplex.Param.mip.pool.relgap.Cur = options.mip.pool.relgap;
    cplex.Param.mip.pool.replace.Cur = options.mip.pool.replace;  
end
temp = cplex.solve();
output.cplexstatus = temp.status;
fval = cplex.Solution.objval;
x = cplex.Solution.x;
if options.mip.pool.intensity
    for i = 2:1:length(temp.pool.solution)
        x = [x temp.pool.solution(i).x];
    end            
end
exitflag = temp.statusstring;
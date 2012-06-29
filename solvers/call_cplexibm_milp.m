function output = call_cplex_milp(interfacedata)

% Author Johan Löfberg
% $Id: call_cplexibm_milp.m,v 1.3 2010-02-08 13:06:11 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
semicont_variables = interfacedata.semicont_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

showprogress(interfacedata.solver.tag,options.showprogress);

if ~isempty(semicont_variables)
    % Bounds must be placed in LB/UB
    [lb,ub,cand_rows_eq,cand_rows_lp] = findulb(F_struc,K,lb,ub);
    F_struc(K.f+cand_rows_lp,:)=[];
    F_struc(cand_rows_eq,:)=[];
    K.l = K.l-length(cand_rows_lp);
    K.f = K.f-length(cand_rows_eq);
    
    redundant = find(lb<=0 & ub>=0);
    semicont_variables = setdiff(semicont_variables,redundant);
end

% Notation used
f = c;
if ~isempty(F_struc)
    Aineq = -F_struc(:,2:end);
    bineq = F_struc(:,1);
else
    Aineq = [];
    bineq = [];
end

if K.f > 0
    Aeq = Aineq(1:K.f,:);
    beq = bineq(1:K.f);
    
    if 0
        Aineq(1:K.f,:)=[];
        bineq(1:K.f)=[];
    else
        [ii,jj,ss] = find(Aineq);keeps = ii>K.f;
        Aineq = sparse(ii(keeps)-K.f,jj(keeps),ss(keeps),size(Aineq,1)-K.f,size(Aineq,2));
        [ii,jj,ss] = find(bineq);keeps = ii>K.f;
        bineq = sparse(ii(keeps)-K.f,jj(keeps),ss(keeps),length(bineq)-K.f,1);
    end   
        
else
    Aeq = [];
    beq = [];
end

% CPLEX assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
NegativeSemiVar = [];
if ~isempty(semicont_variables)
    NegativeSemiVar = find(lb(semicont_variables) < 0);
    if ~isempty(NegativeSemiVar)
        temp = ub(semicont_variables(NegativeSemiVar));
        ub(semicont_variables(NegativeSemiVar)) = -lb(semicont_variables(NegativeSemiVar));
        lb(semicont_variables(NegativeSemiVar)) = -temp;
        if ~isempty(Aineq)
            Aineq(:,semicont_variables(NegativeSemiVar)) = -Aineq(:,semicont_variables(NegativeSemiVar));
        end
        if ~isempty(Aeq)
            Aeq(:,semicont_variables(NegativeSemiVar)) = -Aeq(:,semicont_variables(NegativeSemiVar));
        end
        f(semicont_variables(NegativeSemiVar)) = -f(semicont_variables(NegativeSemiVar));
    end
end

ctype = char(ones(length(f),1)*67);
ctype(binary_variables)  = 'B';  % Should not happen except from bmibnb
if any(semicont_variables)
    ctype(setdiff(integer_variables,semicont_variables)) = 'I';
    ctype(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
    ctype(intersect(semicont_variables,integer_variables)) = 'N';
else
    ctype(integer_variables) = 'I';  
end

% Call mex-interface
if options.savedebug
    save cplexdebug
end
solvertime = clock;

options.cplex.simplex.display = options.verbose;
options.cplex.mip.display = options.verbose;
options.cplex.barrier.display = options.verbose;
if str2num(interfacedata.solver.subversion)>=12.3
    newver = 1;
   if options.verbose
       options.cplex.diagnostics = 'on';
   else
       options.cplex.diagnostics = 'off';
   end
else
    newver = 0;
end

if ~isempty(K.sos.type)
    for i = 1:length(K.sos.weight)
        K.sos.weight{i} = full(K.sos.weight{i}(:));
    end
    if length(K.sos.weight)==1
        K.sos.weight = K.sos.weight{1};
        K.sos.variables = K.sos.variables{1};
    end
end

% Shift the objective (constant term can not be set in CPLEX?)
if isfield(options.cplex,'mip.tolerances.lowercutoff')
    options.cplex.mip.tolerances.lowercutoff = options.cplex.mip.tolerances.lowercutoff-interfacedata.f;
    options.cplex.mip.tolerances.uppercutoff = options.cplex.mip.tolerances.uppercutoff-interfacedata.f;
end
if isempty(integer_variables) && isempty(binary_variables) & isempty(semicont_variables) & isempty(K.sos.type)
    if options.verbose
        [x,fval,exitflag,output,lambda] = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);
    else
        % Some version did not respect silence completely!
        evalc('[x,fval,exitflag,output,lambda] = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);');
    end
else
    if isempty(Aineq) && isempty(Aeq)
        if options.verbose
            [x,fval,exitflag,output] = cplexmilp(f,zeros(1,length(f)),1,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype',x0,options.cplex);
        else
            evalc('[x,fval,exitflag,output] = cplexmilp(f,zeros(1,length(f)),1,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype'',x0,options.cplex); ');
        end
    else
        if options.verbose           
            [x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype',x0,options.cplex);
        else
            evalc('[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype'',x0,options.cplex);');
        end
    end
    lambda = [];
end
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Inconstency in early version of CPLEX
if newver
    the_sign = 1;
else
    the_sign = -1;
end
if ~isempty(lambda)
    D_struc = [the_sign*lambda.eqlin;the_sign*lambda.ineqlin];
else
    D_struc = [];
end

% Gurobi assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
if length(x) == length(f)
    if ~isempty(NegativeSemiVar)
        x(NegativeSemiVar) = -x(NegativeSemiVar);
    end
end

% Check, currently not exhaustive...
switch output.cplexstatus
    case {1,101,102}
        problem = 0;
    case {3,103,106}
        problem = 1; % Infeasible
    case {2,20,21,118}
        problem = 2; % Unbounded
    case {4}
        % not sure if unbounded or infeasible? Let's solve the feasibility
        % problem
        % FIXME: Use conflict refiner...
        options.cplex.Diagnostics = 'off';
        options.cplex.Display = 'off';
        if isempty(integer_variables) & isempty(binary_variables)
            [x,fval,exitflag,output,lambda] = cplexlp(f*0,Aineq,bineq,Aeq,beq,lb,ub,[],options.cplex);
        else
            [x,fval,exitflag,output] = cplexmilp(f*0,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype',[],options.cplex);
            lambda = [];
        end
        if (output.cplexstatus == 3) | (output.cplexstatus == 4)
            problem = 1;
        else
            problem = 2;
        end
    case 119
        problem = 15;
    case {10,11,104,105,107,108,111,112}
        problem = 3; % Iteration/time
    case {5,6,109,110}
        problem = 4; % Numerics
    otherwise
        problem = -1;
end
infostr = yalmiperror(problem,'CPLEX-IBM');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.f = c(:);
    solverinput.Aineq = Aineq;
    solverinput.bineq = bineq;
    solverinput.Aeq = Aeq;
    solverinput.beq = beq;
    solverinput.ctype = ctype;
    solverinput.LB = lb;
    solverinput.UB = ub;
    solverinput.x0 = [];
    solverinput.options = options.cplex;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.x = x;
    solveroutput.fval = fval;
    solveroutput.exitflag = exitflag;
    solveroutput.lambda = lambda;
else
    solveroutput = [];
end

if isempty(x)
    x = zeros(length(f),1);
end

% Standard interface
output.Primal      = x(:);
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;


function [x,fval,exitflag,output,lambda] = cplexlpGATEWAY(f,Aineq,bineq,Aeq,beq,lb,ub,x0,options)

cplex = Cplex('lpex6');
cplex.Model.sense = 'minimize';
bL = [beq;-inf(length(bineq),1)];
bU = [beq;bineq];
if isempty(lb)
    lb = -inf(length(f),1);
end

cplex.addRows(bL, [], bU);
cplex.addCols(f, [Aeq;Aineq], lb,ub);

% Optimize the problem
cplex.solve();

% Write the solution
if ~isfield(cplex.Solution,'x')
    x = zeros(length(lb),1);
    fval = 0;
    lambda.eqlin = [];
    lambda.ineqlin = [];
else
    x = cplex.Solution.x;
    fval = cplex.Solution.objval;
    lambda.eqlin = cplex.Solution.dual(1:length(beq));
    lambda.ineqlin = cplex.Solution.dual(length(beq)+1:end);
end
output.cplexstatus = cplex.Solution.status;
exitflag = output.cplexstatus;
function output = call_cplexibm_miqp(interfacedata)

% Author Johan Löfberg

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

[F_struc,K,lb,ub,semicont_variables] = extractSemiContBounds(F_struc,K,lb,ub,semicont_variables);

% Notation used
H = 2*interfacedata.Q;
f = c;
if ~isempty(F_struc)
    Aineq = -F_struc(:,2:end);
    bineq = F_struc(:,1);
else
    Aineq = [];
    bineq = [];
end

% CPLEX assumes semi-continuous variables only can take negative values so
% we negate semi-continuous violating this
[NegativeSemiVar,H,f,Aineq,lb,ub,semicont_variables] = negateNegativeSemiContVariables(H,f,Aineq,lb,ub,semicont_variables,[]);

if K.f > 0
    Aeq = Aineq(1:K.f,:);
    beq = bineq(1:K.f);
    % Code around performance bugs in older version of MATLAB
    [ii,jj,ss] = find(Aineq);keeps = ii>K.f;
    Aineq = sparse(ii(keeps)-K.f,jj(keeps),ss(keeps),size(Aineq,1)-K.f,size(Aineq,2));
    [ii,jj,ss] = find(bineq);keeps = ii>K.f;
    bineq = sparse(ii(keeps)-K.f,jj(keeps),ss(keeps),length(bineq)-K.f,1);
else
    Aeq = [];
    beq = [];
end

if all(isinf(lb))
    lb = [];
end
if all(isinf(ub))
    ub = [];
end

ctype = char(ones(length(f),1)*67);
ctype(setdiff(integer_variables,semicont_variables)) = 'I';
ctype(binary_variables)  = 'B';  % Should not happen except from bmibnb
ctype(setdiff(semicont_variables,integer_variables)) = 'S';  % Should not happen except from bmibnb
ctype(intersect(semicont_variables,integer_variables)) = 'N';

solvertime = clock;
options.cplex.simplex.display = options.verbose;
options.cplex.mip.display = options.verbose;
options.cplex.barrier.display = options.verbose;
if str2num(interfacedata.solver.subversion)>=12.3
    if options.verbose
        options.cplex.diagnostics = 'on';
    else
        options.cplex.diagnostics = 'off';
    end
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

if options.savedebug
    save cplexdebug
end

% Bug in older versions of CPLEX
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
            evalc('[x,fval,exitflag,output,lambda] = cplexqp(H,f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);');
        else
            evalc('[x,fval,exitflag,output,lambda] = cplexlp(f,Aineq,bineq,Aeq,beq,lb,ub,x0,options.cplex);');
        end
    end
    if ~isempty(lambda) & fixedAineqBug
        lambda.ineqlin = [];
    end
else
    if options.verbose
        if isempty(H)
            [x,fval,exitflag,output] = cplexmilp(f,zeros(1,length(f)),1,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype',x0,options.cplex);
        else
            [x,fval,exitflag,output] = cplexmiqp(H,f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype',x0,options.cplex);
        end
    else
        if isempty(H)
            evalc('[x,fval,exitflag,output] = cplexmilp(f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype'',x0,options.cplex);');
            
        else
            evalc('[x,fval,exitflag,output] = cplexmiqp(H,f,Aineq,bineq,Aeq,beq,K.sos.type,K.sos.variables,K.sos.weight,lb,ub,ctype'',x0,options.cplex);');
        end
    end
    lambda = [];
end
if interfacedata.getsolvertime solvertime = etime(clock,solvertime);else solvertime = 0;end

% Inconstency in early version of CPLEX
if str2num(interfacedata.solver.subversion)>=12.3
    the_sign = 1;
else
    the_sign = -1;
end
if ~isempty(lambda)
    D_struc = [the_sign*lambda.eqlin;the_sign*lambda.ineqlin];
else
    D_struc = [];
    lambda = [];
end

if length(x) == length(f)
    if ~isempty(NegativeSemiVar)
        x(NegativeSemiVar) = -x(NegativeSemiVar);
    end
end

showprogress('Calling CPLEX-IBM',options.showprogress);
% Check, currently not exhaustive...
switch output.cplexstatus
    case {1,101,102}
        problem = 0;
    case {3,103,106}
        problem = 1; % Infeasible
    case {2,20,21,118}
        problem = 2; % Unbounded
    case 4
        if isempty(H)
            options.cplex.Diagnostics = 'off';
            options.cplex.Display = 'off';
            if isempty(integer_variables) & isempty(binary_variables)
                [x,fval,exitflag,output,lambda] = cplexlp(f*0,Aineq,bineq,Aeq,beq,lb,ub,[],options.cplex);
            else
                [x,fval,exitflag,output] = cplexmilp(f*0,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype',[],options.cplex);
                lambda = [];
            end            
        else
            if isempty(integer_variables) & isempty(binary_variables)
                if isempty(Aineq) & isempty(Aeq) & isempty(lb) & isempty(ub)
                    [x,fval,exitflag,output,lambda] = cplexqp(H*0,f*0,zeros(1,length(f)),1,Aeq,beq,lb,ub,[],options.cplex);
                    if ~isempty(lambda)
                        lambda.ineqlin = [];
                    end
                else
                    [x,fval,exitflag,output,lambda] = cplexqp(H*0,f*0,Aineq,bineq,Aeq,beq,lb,ub,[],options.cplex);
                end
            else
                if isempty(Aineq) & isempty(Aeq) & isempty(lb) & isempty(ub)
                    [x,fval,exitflag,output] = cplexmiqp(H*0,f*0,zeros(1,length(f)),1,Aeq,beq,[],[],[],lb,ub,ctype',[],options.cplex);
                else
                    [x,fval,exitflag,output] = cplexmiqp(H*0,f*0,Aineq,bineq,Aeq,beq,[],[],[],lb,ub,ctype',[],options.cplex);
                end
                lambda = [];
            end
        end
        if output.cplexstatus == 3
            problem = 1;
        else
            problem = 2;
        end       
    case {10,11,104,105,107,108,111,112}
        problem = 3; % Iteration/time
    case 119
        problem = 15;
    case {5,6,109,110}
        problem = 4; % Numerics
    otherwise
        problem = -1;
end

infostr = yalmiperror(problem,'CPLEX-IBM');

% Save all data sent to solver?
if options.savesolverinput
    solverinput.H = H;
    solverinput.f = f;
    solverinput.Aineq = Aineq;
    solverinput.Aineq = Aeq;
    solverinput.bineq = bineq;
    solverinput.beq = beq;
    solverinput.LB = lb;
    solverinput.UB = ub;
    solverinput.X0 = [];
    solverinput.ctype = ctype;
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
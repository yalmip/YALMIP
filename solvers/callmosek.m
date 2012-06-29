function output = callmosek(interfacedata)

% Author Johan Löfberg 
% $Id: callmosek.m,v 1.18 2009-03-11 09:45:32 joloef Exp $

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
Q       = interfacedata.Q;
K       = interfacedata.K;
x0      = interfacedata.x0;
integer_variables = interfacedata.integer_variables;
binary_variables = interfacedata.binary_variables;
extended_variables = interfacedata.extended_variables;
ub      = interfacedata.ub;
lb      = interfacedata.lb;
mt      = interfacedata.monomtable;

% *********************************
% What type of variables do we have
% *********************************
linear_variables = find((sum(abs(mt),2)==1) & (any(mt==1,2)));
nonlinear_variables = setdiff((1:size(mt,1))',linear_variables);
sigmonial_variables = find(any(0>mt,2) | any(mt-fix(mt),2));

if ~isempty(sigmonial_variables) | isequal(interfacedata.solver.version,'GEOMETRIC')
    [x,D_struc,problem,res,solvertime,prob] = call_mosek_geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);          
else
    [x,D_struc,problem,res,solvertime,prob] = call_mosek_lpqp(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables,x0);        
end

infostr = yalmiperror(problem,'MOSEK');	

% Save all data sent to solver?
if options.savesolverinput
    solverinput.prob = prob;
else
    solverinput = [];
end

% Save all data from the solver?
if options.savesolveroutput
    solveroutput.res = res;
else
    solveroutput = [];
end

% Standard interface 
output.Primal      = x;
output.Dual        = D_struc;
output.Slack       = [];
output.problem     = problem;
output.infostr     = infostr;
output.solverinput = solverinput;
output.solveroutput= solveroutput;
output.solvertime  = solvertime;


function [x,D_struc,problem,res,solvertime,prob] = call_mosek_lpqp(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables,x0);

prob.c = c;
if ~isempty(F_struc)
    prob.a = -F_struc(:,2:end);
    prob.buc = full(F_struc(:,1));
    prob.blc = repmat(-inf,length(prob.buc),1);
else
    prob.a = sparse(ones(1,length(c))); % Dummy constraint
    prob.buc = inf;
    prob.blc = -inf;
end
if isempty(lb)
    prob.blx = repmat(-inf,1,length(c));
else
    prob.blx = lb;
end
if isempty(ub)
    prob.bux = repmat(inf,1,length(c));
else
    prob.bux = ub;
end


if K.f>0
    prob.blc(1:K.f) = prob.buc(1:K.f);
end

[prob.qosubi,prob.qosubj,prob.qoval] = find(tril(sparse(2*Q)));    

if K.q(1)>0
    nof_new = sum(K.q);
    prob.a = [prob.a [spalloc(K.f,nof_new,0);
              spalloc(K.l,nof_new,0);speye(nof_new)]];
    prob.blc(1+K.f+K.l:end) = prob.buc(1+K.f+K.l:end);
    prob.buc(1+K.f+K.l:end) = prob.buc(1+K.f+K.l:end);    
    prob.c = [prob.c;zeros(nof_new,1)];
    top = size(F_struc,2)-1;
    for i = 1:length(K.q)
        prob.cones{i}.type = 'MSK_CT_QUAD';
        prob.cones{i}.sub  = top+1:top+K.q(i);
        prob.blx(top+1:top+K.q(i)) = -inf;
        prob.bux(top+1:top+K.q(i)) = inf;        
        top = top + K.q(i);
    end    
end

if ~isempty(integer_variables)
    prob.ints.sub = integer_variables;
end

param = options.mosek;
if isfield(param,'param')
    param = rmfield(param,'param');
end

if ~isempty(x0)
    if options.usex0
        prob.sol.int.xx = zeros(max([length(Q) size(prob.a,2)]),1);
        prob.sol.int.xx(integer_variables) = x0(integer_variables);
        evalc('[r,res] = mosekopt (''symbcon'')');
        sc = res.symbcon ;
        param.MSK_IPAR_MIO_CONSTRUCT_SOL = sc.MSK_ON;
    end
end

% Debug?
if options.savedebug
    save mosekdebug prob param
end

% Call MOSEK
solvertime = clock; 
showprogress('Calling MOSEK',options.showprogress);
if options.verbose == 0
    [r,res] = mosekopt('minimize echo(0)',prob,param); 
else
    [r,res] = mosekopt('minimize',prob,param);     
end
solvertime = etime(clock,solvertime);

if (r == 1010) | (r == 1011) | (r==1001)
    problem = -5;
    x = [];
    D_struc = [];
elseif r == 1200
    problem = 7;
    x = [];
    D_struc = [];
else
    % Recover solutions
    sol = res.sol;
    if isempty(integer_variables)
        x = sol.itr.xx(1:length(c)); % Might have added new ones                
        D_struc = (sol.itr.suc-sol.itr.slc);        
        error_message = sol.itr.prosta;
    else        
        try
        x = sol.int.xx(1:length(c)); % Might have added new ones
        D_struc = [];
        error_message = sol.int.prosta;
        catch
            x = [];
            error_message = 'crash';
            D_struc = [];
        end
      
    end
    
    switch error_message
    case {'PRIMAL_AND_DUAL_FEASIBLE','PRIMAL_FEASIBLE'}
        problem = 0;
    case 'PRIMAL_INFEASIBLE'
        problem = 1;
    case 'DUAL_INFEASIBLE'        
        problem = 2;
    case 'PRIMAL_INFEASIBLE_OR_UNBOUNDED'
        problem = 12;
    otherwise
        problem = -1;
    end
end



function [x,D_struc,problem,res,solvertime,prob] = call_mosek_qclp(options,F_struc,c,Q,K,ub,lb,mt,linear_variables);

nonlinear_variables =setdiff(1:length(c),linear_variables);
var_prod = [];
for i = 1:length(nonlinear_variables)
    vars = find(mt(nonlinear_variables(i),:));
    if length(vars)==1
        var_prod = [var_prod;vars vars];
    else
        var_prod = [var_prod;sort(vars)];
    end
end

A = -F_struc(:,1+linear_variables);
u_c = full(F_struc(:,1));
l_c = repmat(-inf,length(u_c),1);
l_x = lb;
u_x = ub;

if K.f>0
    l_c(1:K.f) = u_c(1:K.f);
end

prob.c = c(linear_variables);;
qcsubk = [];
qcsubi = [];
qcsubj = [];
qcval = [];
for k = 1:size(F_struc,1)
    nonlins = find(F_struc(k,nonlinear_variables+1));
    for i = 1:length(nonlins)
        qcsubk = [qcsubk;k];
        qcsubi = [qcsubi;var_prod(nonlins(i),1)];
        qcsubj = [qcsubj;var_prod(nonlins(i),2)];
        if var_prod(nonlins(i),2)==var_prod(nonlins(i),1)
            qcval = [qcval;-2*F_struc(k,1+nonlinear_variables(nonlins(i)))];
        else
            qcval = [qcval;-F_struc(k,1+nonlinear_variables(nonlins(i)))];
        end
    end
end
prob.qcsubk = qcsubk;
prob.qcsubi = qcsubi;
prob.qcsubj = qcsubj;
prob.qcval = qcval;
prob.a = A;
prob.buc = u_c;

param = options.mosek;
if isfield(param,'param')
    param = rmfield(param,'param');
end

% For debugging
if options.savedebug
    save mosekdebug prob param
end

% Call MOSEK
solvertime = clock; 
showprogress('Calling MOSEK',options.showprogress);
solvertime = clock; 
if options.verbose == 0
    [r,res] = mosekopt('minimize echo(0)',prob,param); 
else
    [r,res] = mosekopt('minimize',prob,param);     
end
solvertime = etime(clock,solvertime);

x = zeros(length(c),1);
x(linear_variables)=res.sol.itr.xx;
D_struc = res.sol.itr.suc;
if K.f>0
    D_struc(1:K.f) = (sol.itr.suc(1:K.f)-sol.itr.slc(1:K.f));
end
% Check, currently not exhaustive...
switch res.sol.itr.prosta
case 'PRIMAL_AND_DUAL_FEASIBLE'
    problem = 0;
case 'PRIMAL_INFEASIBLE'
    problem = 1;
case 'DUAL_INFEASIBLE'
    problem = 2;
otherwise
    problem = -1;
end

function [x,D_struc,problem,res,solvertime,prob] = call_mosek_geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);
   
x = [];
D_struc = [];
res = [];
solvertime = 0;
[prob,problem] = yalmip2geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);
if problem == 0

    % Mosek does not support equalities
    if ~isempty(prob.G)
        prob.A = [prob.A;prob.G;-prob.G];
        prob.b = [prob.b;prob.h;1./prob.h];
        prob.map = [prob.map;max(prob.map) + (1:2*length(prob.h))'];
    end
    if options.savedebug
        save mosekdebug prob
    end

    param = options.mosek;
    if isfield(param,'param')
        param = rmfield(param,'param');
    end

    % Call MOSEK
    solvertime = clock;
    showprogress('Calling MOSEK',options.showprogress);
    solvertime = clock;
    if options.verbose == 0  
        res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize echo(0)');
    else
    	res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize');     
    end
    sol = res.sol;
    solvertime = etime(clock,solvertime);

    x = zeros(length(c),1);
    x(linear_variables)=exp(res.sol.itr.xx);
    D_struc = [];%exp(res.sol.itr.suc);

    % Check, currently not exhaustive...
    switch res.sol.itr.prosta
        case 'PRIMAL_AND_DUAL_FEASIBLE'
            problem = 0;
        case 'PRIMAL_INFEASIBLE'
            problem = 1;
        case 'DUAL_INFEASIBLE'
            problem = 2;
        case 'UNKNOWN'
            problem = 9;
        otherwise
            problem = -1;
    end
end
function [x,D_struc,problem,r,res,solvertime,prob] = call_mosek_geometric(model);
   
x = [];
D_struc = [];
r = [];
res = [];
solvertime = 0;
[prob,problem] = yalmip2geometric(model.options,model.F_struc,model.c,model.Q,model.K,model.ub,model.lb,model.monomtable,model.linear_variables,model.extended_variables);
if problem == 0

    % Mosek does not support equalities
    if ~isempty(prob.G)
        prob.A = [prob.A;prob.G;-prob.G];
        prob.b = [prob.b;prob.h;1./prob.h];
        prob.map = [prob.map;max(prob.map) + (1:2*length(prob.h))'];
    end
    if model.options.savedebug
        save mosekdebug prob
    end

    param = model.options.mosek;

    % Call MOSEK   
    showprogress('Calling MOSEK',model.options.showprogress);    
    if model.options.verbose == 0  
        solvertime = tic;
        try
            res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize echo(0)');
        catch
            res = msk9_layer(prob.b,prob.A,prob.map,param,0);
        end
        solvertime = toc(solvertime);
    else
        solvertime = tic;
    	try
            res = mskgpopt(prob.b,prob.A,prob.map,param,'minimize');     
        catch
            res = msk9_layer(prob.b,prob.A,prob.map,param,1);
        end
        solvertime = toc(solvertime);
    end
    sol = res.sol;
    
    x = zeros(length(model.c),1);
    x(model.linear_variables)=exp(res.sol.itr.xx);
    D_struc = [];

    % Check, currently not exhaustive...
    switch res.sol.itr.prosta
        case 'PRIMAL_AND_DUAL_FEASIBLE'
            problem = 0;
        case 'PRIMAL_INFEASIBLE'
            problem = 1;
        case 'DUAL_INFEASIBLE'
            problem = 2;
        case 'UNKNOWN'
            problem = 4;
        otherwise
            problem = -1;
    end
end

function res = msk9_layer(c,A,map,param,verbosity)
% Temporary fix until I do this in low-level format
% Ugly indeed
obj = 0;
con = [];
x = sdpvar(size(A,2),1);
j = find(map == 0);
obj = sum(c(j).*exp(A(j,:)*x));
for i = 1:max(map) 
    j = find(map == i);
    con = [con, sum(c(j).*exp((A(j,:)*x))) <= 1];   
end
sol = optimize(con,obj,sdpsettings('solver','mosek','verbose',verbosity,'savesolveroutput',1));
res = sol.solveroutput.res;
res.sol.itr.xx = value(x);


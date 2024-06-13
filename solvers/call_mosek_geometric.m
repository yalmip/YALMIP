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
        s = 'minimize echo(0)';
    else
        s = 'minimize';
    end
    
    solvertime = tic;
    try
        res = mskgpopt(prob.b,prob.A,prob.map,param,s);
        x = zeros(length(model.c),1);
        x(model.linear_variables)=exp(res.sol.itr.xx);
    catch
        [res,xx] = msk9_layer(prob.b,prob.A,prob.map,param,model);
        x = zeros(length(model.c),1);
        x(model.linear_variables) = xx;
    end
    solvertime = toc(solvertime);
    
    sol = res.sol;
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

function [res,xx] = msk9_layer(c,A,map,param,model)

% First, all variables are re-mapped to x = exp(y) for new variables y
% Data  minimize sum c_j exp(A(j,:)*y) s.t sum c_iK exp(A(K,:))<=1
% We introduce new variables z to represent A*y, and t to represet exp(z)

n_original = size(A,2);
n_new_z = size(A,1);
n_new_t = size(A,1);
model.c = [zeros(n_original,1);zeros(n_new_z,1);c(find(map == 0))];
model.c = [model.c;zeros(n_original+2*n_new_z-length(model.c),1)];

% Create Ay == z
model.F_struc = [spalloc(n_new_z,1,0) A -speye(n_new_z) spalloc(n_new_z,n_new_z,0)];
model.K.f = n_new_z;

% Add all constraints sum c_iK tK <=1
model.K.l = max(map);
for i = 1:max(map)
    k = find(map == i);
    row = sparse(1,k,c(k),1,n_new_t);
    model.F_struc = [ model.F_struc;1 sparse(1,n_original+n_new_z) -row];
end

% Add all expcone definitions exp(z_i) <= t_i
% Format x is EXPCONE <-> x(2)*exp(x(1)/x(2)) <= x(3)
% i.e. [z_i;1;t_i] in EXPCONE
model.K.e = n_new_t;
for i = 1:n_new_t
    rows = sparse(3,n_original+n_new_z*2+1,0);
    rows(1,1 + n_original + i) = 1;
    rows(2,1) = 1;
    rows(3,1 + n_original + n_new_z + i) = 1;
    model.F_struc = [model.F_struc;rows];
end

% Clean model and switch to conic solver
model.variabletype = zeros(1,n_original+n_new_z+n_new_t);
model.monomtable = speye(n_original+n_new_z+n_new_t);
model.sigmonial_variables = [];
model.solver.version = 'SOCP';
model.Q = spalloc(n_original+n_new_z+n_new_t,n_original+n_new_z+n_new_t,0);

% Call recursively with expcone model
model.lb = [];
model.ub = [];
model.options.savesolveroutput = 1;
model.x0 = [];
output = callmosek(model);

% The gluing back to intended problem is a bit messy
res = output.solveroutput.res;
xx = exp(output.Primal(1:n_original));
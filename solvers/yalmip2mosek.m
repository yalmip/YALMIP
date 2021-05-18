function prob = yalmip2mosek(interfacedata);

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
    prob = create_mosek_geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);
else
    prob = create_mosek_lpqp(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables);
end

function prob = create_mosek_lpqp(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,integer_variables);

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

if any(K.q)
    nof_new = sum(K.q);
    prob.a = [prob.a [spalloc(K.f,nof_new,0);spalloc(K.l,nof_new,0);speye(nof_new)]];
    %prob.a(1+K.f+K.l:end,1:length(c)) = prob.a(1+K.f+K.l:end,1:length(c));
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

prob.param = options.mosek;

function prob = create_mosek_geometric(options,F_struc,c,Q,K,ub,lb,mt,linear_variables,extended_variables);

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
    
else
    prob = [];
end
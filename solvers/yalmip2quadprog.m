function model = yalmip2quadprog(interfacedata);

options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

if nnz(Q)==0
    ops = options.linprog;
else
    ops = options.quadprog;
end

switch options.verbose
case 0
    ops.Display = 'off';
case 1
    ops.Display = 'final';
otherwise
    ops.Display = 'iter';
end
    
% QUAPROG does not like lb==ub in (LINUX/6.1)
equality_in_bound = find((abs(lb-ub)<1e-12) & ~isinf(lb));
n = length(c);
m = length(equality_in_bound);
if ~isempty(equality_in_bound)
    F_struc = [-lb(equality_in_bound) sparse(1:m,equality_in_bound,ones(m,1),m,n);F_struc];
    ub(equality_in_bound) = ub(equality_in_bound) + 1;
    lb(equality_in_bound) = lb(equality_in_bound) - 1;
    K.f = K.f + m;
end

if ~isempty(F_struc)
    Aeq = -F_struc(1:1:K.f,2:end);
    beq = F_struc(1:1:K.f,1);        
    A =-F_struc(K.f+1:end,2:end);
    b = F_struc(K.f+1:end,1);   
else
    A = [];
    b = [];
    Aeq = [];
    beq = [];
end

if isfield(ops,'LargeScale')
    if ~isequal(ops.LargeScale,'on')
        Q = full(Q);
        c = full(c);
        A = full(A);
        b = full(b);
        Aeq = full(Aeq);
        beq = full(beq);
    end
end

model.Q = 2*Q;
model.c = c;
model.A = A;
model.b = b;
model.Aeq = Aeq;
model.beq = beq;
model.lb = lb;
model.ub = ub;
model.ops = ops;
model.x0 = interfacedata.x0;
model.f = interfacedata.f;
function model = yalmip2intlinprog(interfacedata)

options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

ops = options.intlinprog;
switch options.verbose
    case 0
        ops.Display = 'off';
    case 1
        ops.Display = 'final';
    otherwise
        ops.Display = 'iter';
end

if isempty(F_struc)
    Aeq = [];
    beq = [];
    A = [];
    b = [];
else
    Aeq = -F_struc(1:1:K.f,2:end);
    beq = F_struc(1:1:K.f,1);        
    A =-F_struc(K.f+1:end,2:end);
    b = F_struc(K.f+1:end,1);   
end

intcon = union(interfacedata.integer_variables,interfacedata.binary_variables);
if ~isempty(interfacedata.binary_variables)
    lb(interfacedata.binary_variables) = max(0,lb(interfacedata.binary_variables));
    ub(interfacedata.binary_variables) = min(1,ub(interfacedata.binary_variables));
end

model.c = c;
model.intcon = intcon;
model.A = A;
model.b = b;
model.Aeq = Aeq;
model.beq = beq;
model.lb = lb;
model.ub = ub;
model.ops = ops;
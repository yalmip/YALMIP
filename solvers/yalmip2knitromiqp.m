function model = yalmip2knitromiqp(interfacedata);

options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

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

ops = options.knitro;
switch options.verbose
    case 0
        ops.Display = 'off';
    otherwise
        ops.Display = 'iter';
end
if isempty(ops.OutputFcn)
    ops=rmfield(ops,'OutputFcn');
end
if isempty(ops.UseParallel)
    ops=rmfield(ops,'UseParallel');
end

xType = zeros(length(interfacedata.c),1);
xType(interfacedata.binary_variables) = 2;
xType(interfacedata.integer_variables) = 1;

model.xType = xType;
model.H = 2*Q;
model.f = c;
model.A = A;
model.b = full(b);
model.Aeq = Aeq;
model.beq = full(beq);
model.lb = lb;
model.ub = ub;
model.ops = ops;
model.x0 = x0;
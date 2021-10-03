function model = yalmip2cbc(interfacedata)

F_struc = interfacedata.F_struc;
Q       = interfacedata.Q;
c       = interfacedata.c;
K       = interfacedata.K;
lb      = full(interfacedata.lb);
ub      = full(interfacedata.ub);
x0      = interfacedata.x0;

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

ivars = repmat('C',length(c),1);
ivars(interfacedata.integer_variables) = 'I';
ivars(interfacedata.binary_variables) = 'B';

% CBC merges equalities and equalities
ru = full([beq;b]);
rl = full([beq;repmat(-inf,length(b),1)]);
A = [Aeq;A];

% SOS currently not supported
sos.type='';
sos.index=[];
sos.weight=[];

model.f = full(c);
model.A = A;
model.rl = rl;
model.ru = ru;
model.lb = lb;
model.ub = ub;
model.xtype = ivars;
model.sos = sos;
model.x0 = [];
model.H = [];
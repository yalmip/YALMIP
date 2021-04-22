function model = yalmip2coneprog(interfacedata);

options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
Q       = interfacedata.Q;
lb      = interfacedata.lb;
ub      = interfacedata.ub;

ops = options.coneprog;
switch options.verbose
case 0
    ops.Display = 'off';
case 1
    ops.Display = 'final';
otherwise
    ops.Display = 'iter';
end
    
Aeq = -F_struc(1:1:K.f,2:end);
beq = F_struc(1:1:K.f,1);        
A =-F_struc(K.f+1:K.f+K.l,2:end);
b = F_struc(K.f+1:K.f+K.l,1);   

if any(K.q)
    top = K.f + K.l + 1;
    for i = 1:length(K.q)
        base = F_struc(top:top+K.q(i)-1,:);
        top = top + K.q(i);
        socConstraints(i) = secondordercone(base(2:end,2:end),-base(2:end,1),base(1,2:end),-base(1));
    end        
else
    socConstraints = [];
end

model.f = c;
model.A = A;
model.b = b;
model.Aeq = Aeq;
model.beq = beq;
model.lb = lb;
model.ub = ub;
model.ops = ops;
model.x0 = interfacedata.x0;
model.socConstraints = socConstraints;
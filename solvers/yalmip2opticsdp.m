function model = sedumi2opticsdp(interfacedata)

% Retrieve needed data
options = interfacedata.options;
F_struc = interfacedata.F_struc;
c       = interfacedata.c;
K       = interfacedata.K;
x0      = interfacedata.x0;
ub      = interfacedata.ub;
lb      = interfacedata.lb;

if interfacedata.K.l > 0
    model.A = -F_struc(1:K.l,2:end);
    model.b = full(F_struc(1:K.l,1));
else
    model.A = [];
    model.b = [];
end

model.lb = lb;
model.ub = ub;
model.f = c;

if K.s(1)>0
    top = 1 + K.l + K.f;
    for j = 1:length(K.s)
        n = K.s(j);
        i = find(triu(ones(n)));
        CA = F_struc(top:top+n^2-1,:);        
        model.sdcone{j} = [-CA(:,1) CA(:,2:end)];
        top = top  + n^2;
    end
else
    model.sdcone = [];
end

model.y0 = x0;
model.ops = options.csdp;

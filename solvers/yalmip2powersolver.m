function model = yalmip2powersolver(interfacedata);

model.A = -interfacedata.F_struc(:,2:end);
model.b = -interfacedata.c;
model.C = interfacedata.F_struc(:,1);
model.K = interfacedata.K;
pars = interfacedata.options.powersolver;


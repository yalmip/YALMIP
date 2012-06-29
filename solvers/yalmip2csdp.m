function model = yalmip2csdp(interfacedata);

pars = interfacedata.options.csdp;
pars.printlevel=interfacedata.options.verbose;
model.At = -interfacedata.F_struc(:,2:end);
model.b = -interfacedata.c;
model.C = interfacedata.F_struc(:,1);
model.K = interfacedata.K;
model.pars = pars;
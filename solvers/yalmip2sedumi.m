function model = yalmip2sedumi(interfacedata);

model.A = -interfacedata.F_struc(:,2:end);
model.b = -interfacedata.c;
model.C = interfacedata.F_struc(:,1);
model.K = interfacedata.K;
pars = interfacedata.options.sedumi;
pars.fid = double(interfacedata.options.verbose);
model.pars = pars;


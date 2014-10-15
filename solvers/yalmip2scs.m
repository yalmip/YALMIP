function model = yalmip2scs(interfacedata);

model.data.A = -interfacedata.F_struc(:,2:end);
model.data.b = full(interfacedata.F_struc(:,1));
model.data.c =  interfacedata.c;
model.cones = interfacedata.K;
model.param = interfacedata.options.scs;


function model = yalmip2dsdp(interfacedata);

[C,A,b,blk] = sedumi2dsdp(interfacedata.F_struc,interfacedata.c,interfacedata.K);
interfacedata.options.dsdp.dual_quadratic=spalloc(length(interfacedata.c),length(interfacedata.c),0);
interfacedata.options.dsdp.printyes = (interfacedata.options.verbose>0);
model.A = A;
model.C = C;
model.b = b;
model.options = interfacedata.options.dsdp

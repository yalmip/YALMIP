function model = yalmip2sdpt3(interfacedata)
[blk,A,C,b,oldKs]=sedumi2sdpt3(interfacedata.F_struc(:,1),-interfacedata.F_struc(:,2:end),-interfacedata.c,interfacedata.K,interfacedata.options.sdpt3.smallblkdim);
interfacedata.options.sdpt3.printyes=double(interfacedata.options.verbose);
interfacedata.options.sdpt3.expon=interfacedata.options.sdpt3.expon(1);

model.blk = blk;
model.A = A;
model.C = C;
model.b = b;
model.ops = interfacedata.options.sdpt3;
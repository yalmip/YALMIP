function model = yalmip2sdpa(interfacedata);
[mDIM,nBLOCK,bLOCKsTRUCT,c,F] = sedumi2sdpa(interfacedata.F_struc,interfacedata.c,interfacedata.K);
if interfacedata.options.verbose==0
    interfacedata.options.sdpa.print = 'no';
else
    interfacedata.options.sdpa.print = 'display';
end
model.mDIM = mDIM;
model.nBLOCK = nBLOCK;
model.bLOCKsTRUCT = bLOCKsTRUCT;
model.c = c;
model.F = F;
model.x0 = [];
model.X0 = [];
model.Y0 = [];
model.OPTIONS = interfacedata.options.sdpa;
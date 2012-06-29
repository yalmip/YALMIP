function  model = yalmip2glpkmex(interfacedata);

n = length(interfacedata.c);
LB=repmat(-1e6,n,1);
UB=repmat(1e6,n,1);
SENSE = 1;
C = full(interfacedata.c);
A =-interfacedata.F_struc(:,2:end);
B = full(interfacedata.F_struc(:,1));
if length(B)==0;
    A = C';
    B = 1e6;
end
CTYPE = [repmat('S',interfacedata.K.f,1); repmat('U',interfacedata.K.l,1)];
VARTYPE = repmat('C',n,1);
VARTYPE(interfacedata.integer_variables)='I';
interfacedata.options.glpk.msglev = interfacedata.options.verbose;
if interfacedata.options.glpk.msglev==1
    interfacedata.options.glpk.msglev = 2;
end
model.SENSE = SENSE;
model.C= C;
model.A = A;
model.B = B;
model.CTYPE = CTYPE;
model.LB = LB;
model.UB = UB;
model.VARTYPE = VARTYPE;
model.PARAM = interfacedata.options.glpk.lpsolver;
model.SAVE = interfacedata.options.glpk.save;
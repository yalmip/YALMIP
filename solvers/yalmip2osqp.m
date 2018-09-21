function model = yalmip2osqp(interfacedata);

% Construct equalities and inequalities
K = interfacedata.K;
F_struc = interfacedata.F_struc;
if ~isempty(F_struc)
    Aeq = -F_struc(1:1:K.f,2:end);
    beq = F_struc(1:1:K.f,1);
    A =-F_struc(K.f+1:end,2:end);
    b = F_struc(K.f+1:end,1);
else
    A = [];
    b = [];
    Aeq = [];
    beq = [];
end

model.options = interfacedata.options.osqp;
model.P = 2*interfacedata.Q;
model.q = interfacedata.c;
eye_n = speye(length(model.q));
model.A = [Aeq; A; eye_n];
model.l = full([beq; -inf(length(b),1); interfacedata.lb]);
model.u = full([beq; b; interfacedata.ub]);

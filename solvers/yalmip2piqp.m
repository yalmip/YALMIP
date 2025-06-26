function model = yalmip2piqp(interfacedata);

quadprog_model = yalmip2quadprog(interfacedata);

model.options = interfacedata.options.piqp;
model.P = quadprog_model.Q;
model.c = quadprog_model.c;
model.A = quadprog_model.Aeq;
model.b = quadprog_model.beq;
model.G = quadprog_model.A;
model.h = quadprog_model.b;
model.x_lb = quadprog_model.lb;
model.x_ub = quadprog_model.ub;
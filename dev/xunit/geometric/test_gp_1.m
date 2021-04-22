function tests = test_gp_1
tests = functiontests(localfunctions);

function test1(dummy)

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
t = [t1 t2 t3];
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);
F = ((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1 <= 1);
F = [F, t>=0];
sol = optimize(F,obj);

assert(sol.problem == 0);
assert(all(abs(value(t)-[2 0.5   1.4142]) <= 1e-3));

function test3(dummy)
sdpvar h w d

Awall  = 1;
Afloor = 1;

F = (0.5 <= h/w <= 2) + (0.5 <= d/w <= 2);
F = F + (2*(h*w+h*d) <= Awall) + (w*d <= Afloor);
F = [F, [h w d] >=0];
sol = optimize(F,-(h*w*d))

assert(sol.problem == 0);
assert(abs(value(h*w*d) - 0.19245008957514) <= 1e-5);
assert(all(abs(value([h w d]) - [ 0.28867519677435   0.57735039354827   1.15470006291485]) <= 1e-2));

function test4(dummy)

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);

F = (max((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1,0.25*t1*t2) <= min(t1,t2));
F = [F, [t1 t2 t3] >= 0];
sol = optimize(F,obj);

assert(sol.problem == 0);
assert(all(abs(value([t1 t2 t3]) - [ 1.10978618937192   1.10978618937162   1.57815225707513]) <=  1e-4));
assert(abs(value(obj) - 1.344555694227871e+002) <= 1e-3);

function test5(dummy)

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);

F = (max((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1,0.25*t1*t2) <= min((t1+0.5*t2)^-1,t2));
F = F + ((2*t1+3*t2^-1)^0.5 <= 2);
F = [F, [t1 t2 t3] >= 0];

sol = optimize(F,obj);

assert(sol.problem == 0);
assert(abs(value(obj) - 2.359439050512407e+002) <= 1e-3);
assert(all(abs(value([t1 t2 t3]) - [0.76467168678701   1.23304260692267   4.24155022707061]) <=  1e-3));

function test6(dummy)

q = sdpvar(1,1);
F = (q >= 0);
obj = (1+q)^2.5,
sol = optimize(F,obj,sdpsettings('debug',1,'solver','fmincon-geometric'));

assert(sol.problem == 0);
assert(abs(value(obj) - 1) <= 1e-5);

function test7(dummy)
sdpvar x y;
sol = optimize([y^+.1+x^0.7<=-1,x>=0,y>=0],x+y)
assert(sol.problem == 1);

function test8(dummy)
sdpvar x y;
sol = optimize([y^+.1+x^0.7<=-1;x<=-41,x>=0,y>=0],x+y)
assert(sol.problem == 1);

function test9(dummy)
sdpvar x y;
sol = optimize([y^+.1+x^0.7<=-1;y.^1.1<=-2,x>=0,y>=0],x+y)
assert(sol.problem == 1);
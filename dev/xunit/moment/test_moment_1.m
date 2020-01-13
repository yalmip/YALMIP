function tests = test_moment_1
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar x1 x2 x3
obj = -2*x1+x2-x3;
F = [x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>=0;
                  4-(x1+x2+x3)>=0;
                   6-(3*x2+x3)>=0;
         2>=x1>=0,x2>=0,3>=x3>=0];
sol = solvemoment(F,obj);
assert(sol.problem == 0);
assert(abs(value(obj)--6) <= 1e-3)

sol = optimize(F,obj,sdpsettings('solver','moment'));
assert(sol.problem == 0);
assert(abs(value(obj)--6) <= 1e-3)

[sol,x,momentdata] = solvemoment(F,obj,[],4);
assert(abs(value(obj)--4) <= 1e-3)
assert(norm(x{1}-[0.5;0;3]) <= 1e-2)
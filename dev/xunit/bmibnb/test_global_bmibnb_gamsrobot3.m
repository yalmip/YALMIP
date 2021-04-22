function tests = test_global_bmibnb_gamsrobot3
tests = functiontests(localfunctions);

function test1(dummy)

Q = 100*eye(10);
c = [48, 42, 48, 45, 44, 41, 47, 42, 45, 46]';
b = [-4, 22,-6,-23,-12]';
A =[-2 -6 -1 0 -3 -3 -2 -6 -2 -2;
6 -5 8 -3 0 1 3 8 9 -3;
-5 6 5 3 8 -8 9 2 0 -9;
9 5 0 -9 1 -8 3 -9 -9 -3;
-8 7 -4 -5 -9 1 -7 -1 3 -2];
x = sdpvar(10,1);
t = sdpvar(1,1);
p = c'*x-0.5*x'*Q*x;
F = (0<=x<=1)+(A*x<=b);
obj = p;

sol = optimize(F,obj,sdpsettings('solver','bmibnb'))

assert(abs(value(obj)--39) <= 1e-3)
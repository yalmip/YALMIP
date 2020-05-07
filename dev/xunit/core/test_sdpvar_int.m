function tests = test_sdpvar_int
tests = functiontests(localfunctions);

function test1(dummy)

% Ordering should not make difference
sdpvar x1 x2
p = 4*x1^4 + x1*x2;
I1 = int(p,[x2 x1],[0 2],[1 5]);
I2 = int(p,[x1 x2],[2 0],[5 1]);
assert(I1 == I2)

sdpvar x1 x2 u
p = 4*x1^4 + x1*x2;
I1 = int(p,[x2 x1],[0 2],[1 u]);
I2 = int(p,[x1 x2],[2 0],[u 1]);
assert((I1-I2) == 0);

I1 = int(p,[x2 x1]);
I2 = 0.2500*x1^2*x2^2+0.8000*x1^5*x2;
assert((I1-I2) == 0);
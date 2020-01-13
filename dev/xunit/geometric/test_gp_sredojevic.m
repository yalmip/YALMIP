function tests = test_gp_sredojevic
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar x y
css=(x>=1)+(y>=1)+(x/y<=4)+(y<=8);
css=css+(x^2/y==1.5);
obj = x+y/x;
sol = optimize(css,obj)

assert(sol.problem == 0);
assert(abs(value(obj)-2.04124145231932)<=1e-4);
assert(all(abs(value([x y])-[ 1.22474487139159   1.00000000000000])<=1e-4));

function test2(dummy)
sdpvar x y z
css=(x>=1)+(y>=1)+(x/y<=4)+(y<=8) + (x*z == 10) %+ (1<z<16);
css=css+(x^2/y==1.5);
obj = x+y/x;
sol = optimize(css,obj)

assert(sol.problem == 0);
assert(abs(value(obj)-2.04124145231932)<=1e-4);
assert(all(abs(value([x y])-[ 1.22474487139159   1.00000000000000])<=1e-4));


function tests = test_operator_entropy
tests = functiontests(localfunctions);

function test1(dummy)
x = sdpvar(2,1);
assign(x,1);
sol = optimize((x >= 0.1),-entropy(x),sdpsettings('usex0',1,'solver','fmincon'));

assert(sol.problem == 0);
assert(abs(value(sum(x))-0.73) <= 1e-2);
assert(abs(-entropy(value(x))--0.7357588) <= 1e-3);

function test2(dummy)
x1 = sdpvar(1,1);
y = sdpvar(1,1);
x2 = sdpvar(1,1);
x = [x1;x2];

assign(x,1);
sol = optimize(([x;y] >= 0.1),-entropy(x)+y,sdpsettings('usex0',1,'solver','fmincon'));

assert(sol.problem == 0);
assert(abs(value(sum(x))-0.73) <= 1e-2);


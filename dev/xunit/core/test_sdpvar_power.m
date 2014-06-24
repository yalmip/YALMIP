function test_sdpvar_power

x = sdpvar(1);
y = ([1 1 1]*x).^(0:2);
assertEqual(getbase(y), getbase([1 x x^2]));
x = sdpvar(1);
y = x.^(0:2);
assertEqual(getbase(y), getbase([1 x x^2]));


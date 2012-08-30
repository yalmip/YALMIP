function test_kktqp

Q = magic(5);
x = sdpvar(5,1);
solvesdp([-1 <= x <= 1],x'*Q*x,sdpsettings('solver','kktqp'))

mbg_asserttolequal(double(OO), -80.9412, 1e-2);

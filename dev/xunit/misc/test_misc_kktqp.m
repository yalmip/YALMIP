function test_misc_kktqp

Q = magic(5);
x = sdpvar(5,1);
OO = x'*Q*x;
solvesdp([-1 <= x <= 1],x'*Q*x,sdpsettings('solver','kktqp'))
assertElementsAlmostEqual(double(OO),-80.9412,'absolute',1e-2);

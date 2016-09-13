function test_operator_max

% Test for bug #345
a=sdpvar(1,1);
b=sdpvar(1,1);
C=[];
C=[C,a==-1];
C=[C,b==max(a,0)];
sol = optimize(C)
assertElementsAlmostEqual(value(b),0, 'absolute',1e-4);


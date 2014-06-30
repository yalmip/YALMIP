function test_moment_1

sdpvar x1 x2 x3
obj = -2*x1+x2-x3;
F = [x1*(4*x1-4*x2+4*x3-20)+x2*(2*x2-2*x3+9)+x3*(2*x3-13)+24>=0;
                  4-(x1+x2+x3)>=0;
                   6-(3*x2+x3)>=0;
         2>=x1>=0,x2>=0,3>=x3>=0];
sol = solvemoment(F,obj);
assertTrue(sol.problem == 0);
assertElementsAlmostEqual(double(obj),-6,'absolute',1e-3);

sol = solvesdp(F,obj,sdpsettings('solver','moment'));
assertTrue(sol.problem == 0);
assertElementsAlmostEqual(double(obj),-6,'absolute',1e-3);

[sol,x,momentdata] = solvemoment(F,obj,[],4);
assertElementsAlmostEqual(double(obj),-4,'absolute',1e-3);
assertElementsAlmostEqual(x{1},[0.5;0;3],'absolute',1e-3);
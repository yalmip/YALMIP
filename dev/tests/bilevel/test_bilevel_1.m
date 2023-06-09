function tests = test_bilevel
tests = functiontests(localfunctions);

function test1(testCase)
ops = sdpsettings('verbose',0);
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2;
CO = [x1 x2]>=0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1];
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],ops);
testCase.assertTrue(abs(value(OO) - (-27.6)) <= 1e-5);
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',0,'verbose',0));
testCase.assertTrue(abs(value(OO) - (-27.6)) <= 1e-5);

function test2(testCase)
% No y3 in outer
ops = sdpsettings('verbose',0);
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2+x1^2+x2^2;
CO = [x1 x2]>=0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1];
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('quadprog.Algorithm','interior-point-convex','verbose',0));
testCase.assertTrue(abs(value(OO) -(-26.7900)) <= 1e-5);

function test3(testCase)
% Lp from bard, with added quadratic term handled by SDP
sdpvar x1 x2 y1 y2 y3 t1 t2
OO = -8*x1-4*x2+4*y1-40*y2+t1+t2;
CO = [[x1 x2]>=0,[t1 x1;x1 1]>=0,[t2 x2;x2 1]>=0];
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1];
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',1,'verbose',0));
testCase.assertTrue(abs(value(OO) -(-26.7900)) <= 1e-5);
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',0,'verbose',0));
testCase.assertTrue(abs(value(OO) -(-26.7900)) <= 1e-5);

function test4(testCase)
% Lp from Bard, with integer variables
ops = sdpsettings('verbose',0);
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2-4*y3;
CO = [[x1 x2]>=0, integer([x1 x2])];
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 10, 2*x1-y1+2*y2-0.5*y3 <= 10, 2*x2+2*y1-y2-0.5*y3 <= 9.7];
sol = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],ops);
testCase.assertTrue(abs(value(OO) -( -256.2667)) <= 1e-4);

function test5(testCase)
% Audet
ops = sdpsettings('verbose',0);
sdpvar x y 
OO = x+2*y;
CO = x>=0;
OI = -4*x-y;
CI = [x+y>=8,3*x-2*y>=-6, -3*x-4*y>=-48,y>=0];
solvebilevel(CO,OO,CI,OI,y,ops);
testCase.assertTrue(abs(value(OO)-14) <= 1e-5);

function test6(testCase)
ops = sdpsettings('verbose',0);
sdpvar x y
OO = -x-10*y;
CO = [x>=0,y>=0];
OI = y;
CI = [-25*x+20*y<=30,x+2*y<=10,2*x-y <= 15,2*x+10*y>=15,x>=0,y>=0];
[sol,info] = solvebilevel(CO,OO,CI,OI,y,ops);
testCase.assertTrue(abs(value(OO)-(-18)) <= 1e-5);

function test7(testCase)
% Crashed due to bug in compress_evaluation_scheme
% Also crashes if evaluationschem is compressed in solvebilevel
sdpvar x1 x2 y1 y2 y3
OO = -8*x1 - 4*x2 + 4*y1 - 40*y2 + 4*y3;
CO = [x1 x2]>=0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1];
solvebilevel(CO,OO+exp(x1)+exp(x2),CI,OI,[y1 y2 y3],sdpsettings('bilevel.outersolver','fmincon','verbose',0));
testCase.assertTrue(abs(value(OO+exp(x1)+exp(x2)) -(-22.5404)) <= 1e-3);

function test8(testCase)
sdpvar y1 y2 y3
sdpvar x1 x2 
OO = -8*x1 - 4*x2 + 4*y1 - 40*y2 + 4*y3;
CO = [[y1 y2 y3]>=0,[x1 x2]>=0];
OI = x1 + 2*x2 + y1 + y2 + 2*y3;
CI = [-y1 + y2 + y3 <= 1,
      2*x1 - y1 + 2*y2 - 0.5*y3 <= 1,
      2*x2 + 2*y1 - y2 - 0.5*y3 <= 1,
                     [y1 y2 y3] >= 0,
                        [x1 x2] >= 0];
solvebilevel(CO,OO-OO^2,CI,OI,[y1 y2 y3],sdpsettings('bilevel.outersolver','bmibnb','bmibnb.uppersolver','fmincon','bilevel.solvefrp',0,'verbose',0));
testCase.assertTrue(abs(value(OO-OO^2) -(-702)) <= 1e-5);

function test9(testCase)
ops = sdpsettings('verbose',0);
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2-4*y3;
CO = [x1>=0, x2>=0];
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,
       -y1+y2+y3 <= 10,
      2*x1-y1+2*y2-0.5*y3 <= 10,
      2*x2+2*y1-y2-0.5*y3 <= 9.7];
[K,details] = kkt(CI,OI^2,[y1 y2 y3],ops);
ops = sdpsettings('solver','bnb','bnb.solver','quadprog','quadprog.Algorithm','interior-point-convex','verbose',0);
optimize([K,CO,-100<=[x1 x2 y1 y2 y3 details.dual(:)']<=100],OO+OO^2,ops);

testCase.assertTrue(abs(value(OO+OO^2) - (-.25)) <= 1e-5);

function test_bilevel_1


sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2;
CO = [x1 x2]>=0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1]
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3]);
assertElementsAlmostEqual(double(OO), -27.6,'absolute', 1e-5);
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',0));
assertElementsAlmostEqual(double(OO), -27.6,'absolute', 1e-5);

% No y3 in outer
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2+x1^2+x2^2;
CO = [x1 x2]>=0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1]
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('quadprog.Algorithm','interior-point-convex'));
assertElementsAlmostEqual(double(OO), -26.7900,'absolute', 1e-5);

% Lp from bard, with added quadratic term handled by SDP
sdpvar x1 x2 y1 y2 y3 t1 t2
OO = -8*x1-4*x2+4*y1-40*y2+t1+t2;
CO = [[x1 x2]>=0,[t1 x1;x1 1]>=0,[t2 x2;x2 1]>=0];
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1]
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',1));
assertElementsAlmostEqual(double(OO), -26.7900,'absolute', 1e-5);
[sol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',0));
assertElementsAlmostEqual(double(OO), -26.7900,'absolute', 1e-5);

% Lp from Bard, with integer variables
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2-4*y3;
CO = [[x1 x2]>=0, integer([x1 x2])];
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 10, 2*x1-y1+2*y2-0.5*y3 <= 10, 2*x2+2*y1-y2-0.5*y3 <= 9.7]
sol = solvebilevel(CO,OO,CI,OI,[y1 y2 y3]);
assertElementsAlmostEqual(double(OO),  -256.2667,'absolute', 1e-4);

% Audet
sdpvar x y 
OO = x+2*y
CO = x>=0
OI = -4*x-y
CI = [x+y>=8,3*x-2*y>=-6, -3*x-4*y>=-48,y>=0]
solvebilevel(CO,OO,CI,OI,[y]);
assertElementsAlmostEqual(double(OO),  14, 'absolute',1e-5);

sdpvar x y
OO = -x-10*y;
CO = [x>=0,y>=0];
OI = y;
CI = [-25*x+20*y<=30,x+2*y<=10,2*x-y <= 15,2*x+10*y>=15,x>=0,y>=0];
[sol,info] = solvebilevel(CO,OO,CI,OI,[y]);
assertElementsAlmostEqual(double(OO),  -18, 'absolute',1e-5);

% Crashed due to bug in compress_evaluation_scheme
% Also crashes if evaluationschem is compressed in solvebilevel
sdpvar x1 x2 y1 y2 y3
OO = -8*x1 - 4*x2 + 4*y1 - 40*y2 + 4*y3;
CO = [x1 x2]>=0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,-y1+y2+y3 <= 1, 2*x1-y1+2*y2-0.5*y3 <= 1, 2*x2+2*y1-y2-0.5*y3 <= 1]
solvebilevel(CO,OO+exp(x1)+exp(x2),CI,OI,[y1 y2 y3]);
assertElementsAlmostEqual(double(OO+exp(x1)+exp(x2)), -22.5404,'absolute', 1e-3);
function test_bilevel_2

% Nonconvex quadratic. Differs between kktp and bmibnb
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
solvebilevel(CO,OO-OO^2,CI,OI,[y1 y2 y3],sdpsettings('bilevel.outersolver','bmibnb','bmibnb.upper','fmincon','bilevel.solvefrp',0));
assertElementsAlmostEqual(double(OO-OO^2),  -702, 'absolute', 1e-5);
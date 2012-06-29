function test_bilevel_3

sdpvar x1 x2 y1 y2 y3

OO = -8*x1-4*x2+4*y1-40*y2-4*y3;
CO = [x1>=0, x2>=0];

OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] >= 0,
       -y1+y2+y3 <= 10,
      2*x1-y1+2*y2-0.5*y3 <= 10,
      2*x2+2*y1-y2-0.5*y3 <= 9.7];

[K,details] = kkt(CI,OI^2,[y1 y2 y3]);
solvesdp([K,CO,-100<=[x1 x2 y1 y2 y3 details.dual(:)']<=100],OO+OO^2,sdpsettings('solver','bnb'))


mbg_asserttolequal(double(OO+OO^2),  -.25, 1e-5);
function tests = test_geometricprogramming
tests = functiontests(localfunctions);

function test1(testCase)

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
t = [t1 t2 t3];
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);
F = ((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1 <= 1);
F = [F, t>=0];
sol = optimize(F,obj,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(all(abs(value(t)-[2 0.5   1.4142]) <= 1e-3));

function test3(testCase)
sdpvar h w d

Awall  = 1;
Afloor = 1;

F = (0.5 <= h/w <= 2) + (0.5 <= d/w <= 2);
F = F + (2*(h*w+h*d) <= Awall) + (w*d <= Afloor);
F = [F, [h w d] >=0];
sol = optimize(F,-(h*w*d),sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(h*w*d) - 0.19245008957514) <= 1e-5);
testCase.assertTrue(all(abs(value([h w d]) - [ 0.28867519677435   0.57735039354827   1.15470006291485]) <= 1e-2));

function test4(testCase)

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);

F = (max((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1,0.25*t1*t2) <= min(t1,t2));
F = [F, [t1 t2 t3] >= 0];
sol = optimize(F,obj,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(all(abs(value([t1 t2 t3]) - [ 1.10978618937192   1.10978618937162   1.57815225707513]) <=  1e-4));
testCase.assertTrue(abs(value(obj) - 1.344555694227871e+002) <= 1e-3);

function test5(testCase)

t1 = sdpvar(1,1);
t2 = sdpvar(1,1);
t3 = sdpvar(1,1);
obj = (40*t1^-1*t2^-0.5*t3^-1)+(20*t1*t3)+(40*t1*t2*t3);

F = (max((1/3)*t1^-2*t2^-2+(4/3)*t2^0.5*t3^-1,0.25*t1*t2) <= min((t1+0.5*t2)^-1,t2));
F = F + ((2*t1+3*t2^-1)^0.5 <= 2);
F = [F, [t1 t2 t3] >= 0];

sol = optimize(F,obj,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 2.359439050512407e+002) <= 1e-3);
testCase.assertTrue(all(abs(value([t1 t2 t3]) - [0.76467168678701   1.23304260692267   4.24155022707061]) <=  1e-3));

function test6(testCase)

q = sdpvar(1,1);
F = (q >= 0);
obj = (1+q)^2.5;
sol = optimize(F,obj,sdpsettings('verbose',0,'debug',1,'solver','fmincon-geometric'));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj) - 1) <= 1e-5);

function test7(testCase)
sdpvar x y;
sol = optimize([y^+.1+x^0.7<=-1,x>=0,y>=0],x+y,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 1);

function test8(testCase)
sdpvar x y;
sol = optimize([y^+.1+x^0.7<=-1;x<=-41,x>=0,y>=0],x+y,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 1);

function test9(testCase)
sdpvar x y;
sol = optimize([y^+.1+x^0.7<=-1;y.^1.1<=-2,x>=0,y>=0],x+y,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 1);


function test10(testCase)
N = 8;
w = sdpvar(N,1);
h = sdpvar(N,1);

% constants
wmin = .1; wmax = 100;
hmin = .1; hmax = 6;
Smin = 1/5; Smax = 5;
sigma_max = 1;
ymax = 10;
E = 1; F = 1;

% objective is the total volume of the beam
% obj = sum of (widths*heights*lengths) over each section
% (recall that the length of each segment is set to be 1)
obj = w'*h; 

% recursive formulation
v = sdpvar(N+1,1); y = sdpvar(N+1,1);
v(N+1,1) = 0; y(N+1,1) = 0;
for i = N:-1:1
  v(i) = 12*(i-1/2)*F/(E*w(i)*h(i)^3) + v(i+1);
  y(i) = 6*(i-1/3)*F/(E*w(i)*h(i)^3)  + v(i+1) + y(i+1);
end

% constraint set
constr = [ ...
  wmin*ones(N,1) <= w; w <= wmax*ones(N,1);
  hmin*ones(N,1) <= h; h <= hmax*ones(N,1);
  Smin*ones(N,1) <= h./w; h./w <= Smax*ones(N,1);
  6*F*[1:N]'./(w.*(h.^2)) <= sigma_max*ones(N,1);
  y(1) <= ymax;
];

% solve GP and compute the optimal volume
sol  = optimize(constr,obj,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj)-42.39654132455499) <= 1e-3);


function test11(testCase)

x=sdpvar(1,1);
y=sdpvar(1,1);
t=x/y;
F=(t>=0.5);
F=F+(y<=2);
F=F+(y>=1);
obj=x^2*y^3;
F=F+(t<=1);
F=F+(t>=1);
F=F+(y^2<=4);
F = [F, x>=0, y>=0];
sol = optimize(F,obj,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj)-1)<=1e-4);
testCase.assertTrue(all(abs(value([x y t]) - [1 1 1]) <= 1e-4));

sol = optimize(F,1/obj,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj)-32) <= 1e-3);
testCase.assertTrue(all(abs(value([x y t]) - [2 2 1]) <= 1e-4));

function test12(testCase)

x=sdpvar(1,1);
y=sdpvar(1,1);
t=x/y;
F=(t>=0.5);
F=F+(y<=2);
F=F+(y>=1);
obj=x^2*y^3;
F=F+(t<=1);
F=F+(t>=1);
F=F+(y^2.5<=4);
F = [F, x>=0];
sol = optimize(F,obj,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj)-1)<=1e-4);
testCase.assertTrue(all(abs(value([x y t]) - [1 1 1]) <= 1e-4));

sol = optimize(F,1/obj,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj)-16) <= 1e-3);
testCase.assertTrue(all(abs(value([x y t]) - [1.74110112659225   1.74110112659225   1.00000000000000]) <=  1e-4));


function test13(testCase)
sdpvar x y
css=(x>=1)+(y>=1)+(x/y<=4)+(y<=8);
css=css+(x^2/y==1.5);
obj = x+y/x;
sol = optimize(css,obj,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj)-2.04124145231932)<=1e-4);
testCase.assertTrue(all(abs(value([x y])-[ 1.22474487139159   1.00000000000000])<=1e-4));

function test14(testCase)
sdpvar x y z
css=(x>=1)+(y>=1)+(x/y<=4)+(y<=8) + (x*z == 10) %+ (1<z<16);
css=css+(x^2/y==1.5);
obj = x+y/x;
sol = optimize(css,obj,sdpsettings('verbose',0));

testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(obj)-2.04124145231932)<=1e-4);
testCase.assertTrue(all(abs(value([x y])-[ 1.22474487139159   1.00000000000000])<=1e-4));


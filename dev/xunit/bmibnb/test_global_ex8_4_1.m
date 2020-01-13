function tests = test_global_ex8_4_1
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from ex8_4_1.gms
% Created 18-Mar-2008 09:31:40 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
x1 = sdpvar(1);
x2 = sdpvar(1);
x3 = sdpvar(1);
x4 = sdpvar(1);
x5 = sdpvar(1);
x6 = sdpvar(1);
x7 = sdpvar(1);
x8 = sdpvar(1);
x9 = sdpvar(1);
x10 = sdpvar(1);
x11 = sdpvar(1);
x12 = sdpvar(1);
x13 = sdpvar(1);
x14 = sdpvar(1);
x15 = sdpvar(1);
x16 = sdpvar(1);
x17 = sdpvar(1);
x18 = sdpvar(1);
x19 = sdpvar(1);
x20 = sdpvar(1);
x21 = sdpvar(1);
x22 = sdpvar(1);

% Define objective function 
objective = -(-(sqr(x1)+sqr(x2-5.9)+sqr(x3-0.9)+sqr(x4-5.4)+sqr(x5-1.8)+sqr(x6-4.4)+sqr(x7-2.6)+sqr(x8-4.6)+sqr(x9-3.3)+sqr(x10-3.5)+sqr(x11-4.4)+sqr(x12-3.7)+sqr(x13-5.2)+sqr(x14-2.8)+sqr(x15-6.1)+sqr(x16-2.8)+sqr(x17-6.5)+sqr(x18-2.4)+sqr(x19-7.4)+sqr(x20-1.5))+0-(0));

% Define constraints 
F = ([]);
F=[F,x22*x1-x2+x21==0];
F=[F,x22*x3-x4+x21==0];
F=[F,x22*x5-x6+x21==0];
F=[F,x22*x7-x8+x21==0];
F=[F,x22*x9-x10+x21==0];
F=[F,x22*x11-x12+x21==0];
F=[F,x22*x13-x14+x21==0];
F=[F,x22*x15-x16+x21==0];
F=[F,x22*x17-x18+x21==0];
F=[F,x22*x19-x20+x21==0];
F=[F,-0.5<=x1<=0.5];
F=[F,5.4<=x2<=6.4];
F=[F,0.4<=x3<=1.4];
F=[F,4.9<=x4<=5.9];
F=[F,1.3<=x5<=2.3];
F=[F,3.9<=x6<=4.9];
F=[F,2.1<=x7<=3.1];
F=[F,4.1<=x8<=5.1];
F=[F,2.8<=x9<=3.8];
F=[F,3<=x10<=4];
F=[F,3.9<=x11<=4.9];
F=[F,3.2<=x12<=4.2];
F=[F,4.7<=x13<=5.7];
F=[F,2.3<=x14<=3.3];
F=[F,5.6<=x15<=6.6];
F=[F,2.3<=x16<=3.3];
F=[F,6<=x17<=7];
F=[F,1.9<=x18<=2.9];
F=[F,6.9<=x19<=7.9];
F=[F,1<=x20<=2];
F=[F,0<=x21<=10];
F=[F,-2<=x22<=2];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)-0.61857) <= 1e-2)
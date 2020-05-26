function tests = test_global_worst
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from worst.gms
% Created 21-Aug-2007 18:27:04 using YALMIP R20070810

% Setup a clean YALMIP environment 
yalmip('clear') 

% Define all variables 
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
x23 = sdpvar(1);
x24 = sdpvar(1);
x25 = sdpvar(1);
x26 = sdpvar(1);
x27 = sdpvar(1);
x28 = sdpvar(1);
x29 = sdpvar(1);
x30 = sdpvar(1);
x31 = sdpvar(1);
x32 = sdpvar(1);
x33 = sdpvar(1);
x34 = sdpvar(1);
x35 = sdpvar(1);

% Define objective function 
objective = -(0-x18-x19-x20-x21-x22+30000*x23-25000*x24+30000*x25+50000*x26-25000*x27-5000*x28-15000*x29-50000*x30-(20682900));

% Define constraints 
errorf = @(x)((1 + erf(x/sqrt(2)))/2);
F = ([]);
F=[F,-95.54*exp(0.09167*x31)+x18==0];
F=[F,-93.27*exp(0.33889*x32)+x19==0];
F=[F,-95.54*exp(0.09167*x31)+x20==0];
F=[F,-93.27*exp(0.33889*x32)+x21==0];
F=[F,-91.03*exp(0.58889*x33)+x22==0];
F=[F,-exp(-0.33889*x32)*(x21*errorf(x2)-95*errorf(x10))+x23==0];
F=[F,-exp(-0.33889*x32)*(x21*errorf(x3)-97*errorf(x11))+x25==0];
F=[F,-exp(-0.58889*x33)*(x22*errorf(x6)-95*errorf(x14))+x24==0];
F=[F,-exp(-0.58889*x33)*(x22*errorf(x7)-97*errorf(x15))+x26==0];
F=[F,-exp(-0.58889*x33)*(x22*errorf(x8)-99*errorf(x16))+x27==0];
F=[F,-exp(-0.33889*x32)*(95*errorf(-x12)-x21*errorf(-x4))+x28==0];
F=[F,-exp(-0.33889*x32)*(97*errorf(-x13)-x21*errorf(-x5))+x29==0];
F=[F,-exp(-0.58889*x33)*(99*errorf(-x17)-x22*errorf(-x9))+x30==0];
F=[F,-1.71779218689115*(log(0.0105263157894737*x21)+0.169445*sqr(x34))/x34+x2==0];
F=[F,-1.71779218689115*(log(0.0103092783505155*x21)+0.169445*sqr(x34))/x34+x3==0];
F=[F,-1.71779218689115*(log(0.0105263157894737*x21)+0.169445*sqr(x34))/x34+x4==0];
F=[F,-1.71779218689115*(log(0.0103092783505155*x21)+0.169445*sqr(x34))/x34+x5==0];
F=[F,-1.30311549893554*(log(0.0105263157894737*x22)+0.294445*sqr(x35))/x35+x6==0];
F=[F,-1.30311549893554*(log(0.0103092783505155*x22)+0.294445*sqr(x35))/x35+x7==0];
F=[F,-1.30311549893554*(log(0.0101010101010101*x22)+0.294445*sqr(x35))/x35+x8==0];
F=[F,-1.30311549893554*(log(0.0101010101010101*x22)+0.294445*sqr(x35))/x35+x9==0];
F=[F,-x2+x10+0.582142594215541*x34==0];
F=[F,-x3+x11+0.582142594215541*x34==0];
F=[F,-x4+x12+0.582142594215541*x34==0];
F=[F,-x5+x13+0.582142594215541*x34==0];
F=[F,-x6+x14+0.767391686168152*x35==0];
F=[F,-x7+x15+0.767391686168152*x35==0];
F=[F,-x8+x16+0.767391686168152*x35==0];
F=[F,-x9+x17+0.767391686168152*x35==0];
F=[F,0.001<=x18];
F=[F,0.001<=x19];
F=[F,0.001<=x20];
F=[F,0.001<=x21];
F=[F,0.001<=x22];
F=[F,0<=x23];
F=[F,0<=x24];
F=[F,0<=x25];
F=[F,0<=x26];
F=[F,0<=x27];
F=[F,0<=x28];
F=[F,0<=x29];
F=[F,0<=x30];
F=[F,0.05245<=x31<=0.0857];
F=[F,0.06175<=x32<=0.095];
F=[F,0.0619<=x33<=0.0939];
F=[F,0.0368<=x34<=0.0768];
F=[F,0.0368<=x35<=0.0768];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1,'bmibnb.lpreduce',1));
assert(sol.problem==0)
assert(abs(value(objective)-2.0762609e+007) <= 1e3)
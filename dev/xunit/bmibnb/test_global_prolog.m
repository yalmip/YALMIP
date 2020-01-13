function tests = test_global_prolog
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from prolog.gms
% Created 07-Jul-2019 11:28:30 using YALMIP R20190425

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

% Define objective function 
objective = -(x5*x2+x6*x2+x7*x4+x8*x4)-0+3712*x9+5000*x10-(0);

% Define constraints 
F = ([]);
F=[F,x5+x6-0.94*x11-0.94*x12-0.94*x13+0.244*x17+0.244*x18+0.244*x19<=0];
F=[F,0.064*x11+0.064*x12+0.064*x13-0.58*x14-0.58*x15-0.58*x16+0.172*x17+0.172*x18+0.172*x19<=0];
F=[F,x7+x8+0.048*x11+0.048*x12+0.048*x13+0.247*x14+0.247*x15+0.247*x16-0.916*x17-0.916*x18-0.916*x19<=0];
F=[F,x11+1.2*x12+0.8*x13+2*x14+1.8*x15+2.4*x16+3*x17+2.7*x18+3.2*x19<=3712];
F=[F,2*x11+1.8*x12+2.2*x13+3*x14+3.5*x15+2.3*x16+3*x17+3.2*x18+2.7*x19<=5000];
F=[F,356.474947137148*x2+53.7083537310174*x4+x5-0.564264890180399*x20<=352];
F=[F,339.983422262764*x2+43.5418249774113*x4+x6-0.405939876920766*x21<=430];
F=[F,106.946746119538*x2+145.018955433089*x4+x7-0.507117039797071*x20<=222];
F=[F,173.929713444361*x2+203.031384299627*x4+x8-0.578889145413521*x21<=292];
F=[F,x5*x2+x7*x4-x20<=0];
F=[F,x6*x2+x8*x4-x21<=0];
F=[F,-3340.8*x9-500*x10+x20<=0];
F=[F,-371.2*x9-4500*x10+x21<=0];
F=[F,0.94*x2-0.064*x3-0.048*x4-x9-2*x10<=0];
F=[F,0.94*x2-0.064*x3-0.048*x4-1.2*x9-1.8*x10<=0];
F=[F,0.94*x2-0.064*x3-0.048*x4-0.8*x9-2.2*x10<=0];
F=[F,0.58*x3-0.247*x4-2*x9-3*x10<=0];
F=[F,0.58*x3-0.247*x4-1.8*x9-3.5*x10<=0];
F=[F,0.58*x3-0.247*x4-2.4*x9-2.3*x10<=0];
F=[F,-0.244*x2-0.172*x3+0.916*x4-3*x9-3*x10<=0];
F=[F,-0.244*x2-0.172*x3+0.916*x4-2.7*x9-3.2*x10<=0];
F=[F,-0.244*x2-0.172*x3+0.916*x4-3.2*x9-2.7*x10<=0];
F=[F,0.2<=x2];
F=[F,0.2<=x3];
F=[F,0.2<=x4];
F=[F,0<=x5];
F=[F,0<=x6];
F=[F,0<=x7];
F=[F,0<=x8];
F=[F,0<=x9];
F=[F,0<=x10];
F=[F,0<=x11];
F=[F,0<=x12];
F=[F,0<=x13];
F=[F,0<=x14];
F=[F,0<=x15];
F=[F,0<=x16];
F=[F,0<=x17];
F=[F,0<=x18];
F=[F,0<=x19];
F=[F,0<=x20];
F=[F,0<=x21];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem == 0)
assert(abs(value(objective)-0) <= 1e-2)
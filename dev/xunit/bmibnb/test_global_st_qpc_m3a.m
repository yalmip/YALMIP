function tests = test_global_st_qpc_m3a
tests = functiontests(localfunctions);

function test1(dummy)
% Model generated from st_qpc-m3a.gms
% Created 06-Aug-2007 09:27:04 using YALMIP R20070725

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

% Define objective function 
objective = -(-(10*x1-6.8*x1*x1-4.6*x1*x2+10*x2-7.9*x1*x3+10*x3-5.1*x1*x4+10*x4-6.9*x1*x5+10*x5-6.8*x1*x6+10*x6-4.6*x1*x7+10*x7-7.9*x1*x8+10*x8-5.1*x1*x9+10*x9-6.9*x1*x10+10*x10-4.6*x2*x1-5.5*x2*x2-5.8*x2*x3-4.5*x2*x4-6*x2*x5-4.6*x2*x6-5.5*x2*x7-5.8*x2*x8-4.5*x2*x9-6*x2*x10-7.9*x3*x1-5.8*x3*x2-13.3*x3*x3-6.7*x3*x4-8.9*x3*x5-7.9*x3*x6-5.8*x3*x7-13.3*x3*x8-6.7*x3*x9-8.9*x3*x10-5.1*x4*x1-4.5*x4*x2-6.7*x4*x3-6.9*x4*x4-5.8*x4*x5-5.1*x4*x6-4.5*x4*x7-6.7*x4*x8-6.9*x4*x9-5.8*x4*x10-6.9*x5*x1-6*x5*x2-8.9*x5*x3-5.8*x5*x4-11.9*x5*x5-6.9*x5*x6-6*x5*x7-8.9*x5*x8-5.8*x5*x9-11.9*x5*x10-6.8*x6*x1-4.6*x6*x2-7.9*x6*x3-5.1*x6*x4-6.9*x6*x5-6.8*x6*x6-4.6*x6*x7-7.9*x6*x8-5.1*x6*x9-6.9*x6*x10-4.6*x7*x1-5.5*x7*x2-5.8*x7*x3-4.5*x7*x4-6*x7*x5-4.6*x7*x6-5.5*x7*x7-5.8*x7*x8-4.5*x7*x9-6*x7*x10-7.9*x8*x1-5.8*x8*x2-13.3*x8*x3-6.7*x8*x4-8.9*x8*x5-7.9*x8*x6-5.8*x8*x7-13.3*x8*x8-6.7*x8*x9-8.9*x8*x10-5.1*x9*x1-4.5*x9*x2-6.7*x9*x3-6.9*x9*x4-5.8*x9*x5-5.1*x9*x6-4.5*x9*x7-6.7*x9*x8-6.9*x9*x9-5.8*x9*x10-6.9*x10*x1-6*x10*x2-8.9*x10*x3-5.8*x10*x4-11.9*x10*x5-6.9*x10*x6-6*x10*x7-8.9*x10*x8-5.8*x10*x9-11.9*x10*x10)+0-(0));

% Define constraints 
F = ([]);
F=[F,20*x1+20*x2+60*x3+60*x4+60*x5+60*x6+5*x7+45*x8+55*x9+65*x10<=600.1];
F=[F,5*x1+7*x2+3*x3+8*x4+13*x5+13*x6+2*x7+14*x8+14*x9+14*x10<=310.5];
F=[F,100*x1+130*x2+50*x3+70*x4+70*x5+70*x6+20*x7+80*x8+80*x9+80*x10<=1800];
F=[F,200*x1+280*x2+100*x3+200*x4+250*x5+280*x6+100*x7+180*x8+200*x9+220*x10<=3850];
F=[F,2*x1+2*x2+4*x3+4*x4+4*x5+4*x6+2*x7+6*x8+6*x9+6*x10<=18.6];
F=[F,4*x1+8*x2+2*x3+6*x4+10*x5+10*x6+5*x7+10*x8+10*x9+10*x10<=198.7];
F=[F,60*x1+110*x2+20*x3+40*x4+60*x5+70*x6+10*x7+40*x8+50*x9+50*x10<=882];
F=[F,150*x1+210*x2+40*x3+70*x4+90*x5+105*x6+60*x7+100*x8+140*x9+180*x10<=4200];
F=[F,80*x1+100*x2+6*x3+16*x4+20*x5+22*x6+20*x8+30*x9+30*x10<=40.25];
F=[F,40*x1+40*x2+12*x3+20*x4+24*x5+28*x6+40*x9+50*x10<=327];
F=[F,0<=x1];
F=[F,0<=x2];
F=[F,0<=x3];
F=[F,0<=x4];
F=[F,0<=x5];
F=[F,0<=x6];
F=[F,0<=x7];
F=[F,0<=x8];
F=[F,0<=x9];
F=[F,0<=x10];

% Solve problem
sol = optimize(F,objective,sdpsettings('solver','bmibnb','allownonconvex',1));
assert(sol.problem==0)
assert(abs(value(objective)--382.695) <= 1)
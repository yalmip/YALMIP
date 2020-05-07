function tests = test_global_bmibnb_qcqp5
tests = functiontests(localfunctions);

function test1(dummy)

x1 = sdpvar(1,1);
x2 = sdpvar(1,1);
x3 = sdpvar(1,1);
x4 = sdpvar(1,1);
x5 = sdpvar(1,1);
x6 = sdpvar(1,1);
x7 = sdpvar(1,1);
x8 = sdpvar(1,1);
x9 = sdpvar(1,1);
x10 = sdpvar(1,1);
x11 = sdpvar(1,1);
t = sdpvar(1,1);
x = [x1;x2;x3;x4;x5;x6;x7;x8;x9;x10;x11];

% ...and 22 linear constraints
e1 =  - x1 - 2*x2 - 3*x3 - 4*x4 - 5*x5 - 6*x6 - 7*x7 - 8*x8 - 9*x9 - 10*x10      - 11*x11 	- 0;
e2 =   - 2*x1 - 3*x2 - 4*x3 - 5*x4 - 6*x5 - 7*x6 - 8*x7 - 9*x8 - 10*x9 - 11*x10      - x11 	- 0;
e3 =   - 3*x1 - 4*x2 - 5*x3 - 6*x4 - 7*x5 - 8*x6 - 9*x7 - 10*x8 - 11*x9 - x10      - 2*x11 	- 0;
e4 =   - 4*x1 - 5*x2 - 6*x3 - 7*x4 - 8*x5 - 9*x6 - 10*x7 - 11*x8 - x9 - 2*x10      - 3*x11 	- 0;
e5 =   - 5*x1 - 6*x2 - 7*x3 - 8*x4 - 9*x5 - 10*x6 - 11*x7 - x8 - 2*x9 - 3*x10      - 4*x11 	- 0;
e6 =   - 6*x1 - 7*x2 - 8*x3 - 9*x4 - 10*x5 - 11*x6 - x7 - 2*x8 - 3*x9 - 4*x10      - 5*x11 	- 0;
e7 =   - 7*x1 - 8*x2 - 9*x3 - 10*x4 - 11*x5 - x6 - 2*x7 - 3*x8 - 4*x9 - 5*x10      - 6*x11 	- 0;
e8 =   - 8*x1 - 9*x2 - 10*x3 - 11*x4 - x5 - 2*x6 - 3*x7 - 4*x8 - 5*x9 - 6*x10      - 7*x11 	- 0;
e9 =   - 9*x1 - 10*x2 - 11*x3 - x4 - 2*x5 - 3*x6 - 4*x7 - 5*x8 - 6*x9 - 7*x10      - 8*x11 	- 0;
e10 =   - 10*x1 - 11*x2 - x3 - 2*x4 - 3*x5 - 4*x6 - 5*x7 - 6*x8 - 7*x9 - 8*x10       - 9*x11 	- 0;
e11 =   - 11*x1 - x2 - 2*x3 - 3*x4 - 4*x5 - 5*x6 - 6*x7 - 7*x8 - 8*x9 - 9*x10       - 10*x11 	- 0;
e12 =     x1 + 2*x2 + 3*x3 + 4*x4 + 5*x5 + 6*x6 + 7*x7 + 8*x8 + 9*x9 + 10*x10       + 11*x11 	- 66;
e13 =     2*x1 + 3*x2 + 4*x3 + 5*x4 + 6*x5 + 7*x6 + 8*x7 + 9*x8 + 10*x9 + 11*x10       + x11 	- 66;
e14 =     3*x1 + 4*x2 + 5*x3 + 6*x4 + 7*x5 + 8*x6 + 9*x7 + 10*x8 + 11*x9 + x10       + 2*x11 	- 66;
e15 =     4*x1 + 5*x2 + 6*x3 + 7*x4 + 8*x5 + 9*x6 + 10*x7 + 11*x8 + x9 + 2*x10       + 3*x11 	- 66;
e16 =     5*x1 + 6*x2 + 7*x3 + 8*x4 + 9*x5 + 10*x6 + 11*x7 + x8 + 2*x9 + 3*x10       + 4*x11 	- 66;
e17 =     6*x1 + 7*x2 + 8*x3 + 9*x4 + 10*x5 + 11*x6 + x7 + 2*x8 + 3*x9 + 4*x10       + 5*x11 	- 66;
e18 =     7*x1 + 8*x2 + 9*x3 + 10*x4 + 11*x5 + x6 + 2*x7 + 3*x8 + 4*x9 + 5*x10       + 6*x11 	- 66;
e19 =     8*x1 + 9*x2 + 10*x3 + 11*x4 + x5 + 2*x6 + 3*x7 + 4*x8 + 5*x9 + 6*x10       + 7*x11 	- 66;
e20 =     9*x1 + 10*x2 + 11*x3 + x4 + 2*x5 + 3*x6 + 4*x7 + 5*x8 + 6*x9 + 7*x10       + 8*x11 	- 66;
e21 =     10*x1 + 11*x2 + x3 + 2*x4 + 3*x5 + 4*x6 + 5*x7 + 6*x8 + 7*x9 + 8*x10       + 9*x11 	- 66;
e22 =     11*x1 + x2 + 2*x3 + 3*x4 + 4*x5 + 5*x6 + 6*x7 + 7*x8 + 8*x9 + 9*x10       + 10*x11 	- 66;

% Define the objective function and the whole program
obj =  (0.5*x1*x2 - x1*x1 + 0.5*x2*x1 - x2*x2 + 0.5*x2*x3 + 0.5*x3*x2 - x3*x3       + 0.5*x3*x4 + 0.5*x4*x3 - x4*x4 + 0.5*x4*x5 + 0.5*x5*x4 - x5*x5 + 0.5*x5      *x6 + 0.5*x6*x5 - x6*x6 + 0.5*x6*x7 + 0.5*x7*x6 - x7*x7 + 0.5*x7*x8 + 0.5      *x8*x7 - x8*x8 + 0.5*x8*x9 + 0.5*x9*x8 - x9*x9 + 0.5*x9*x10 + 0.5*x10*x9       - x10*x10 + 0.5*x10*x11 + 0.5*x11*x10 - x11*x11);
e = -[e1;e2;e3;e4;e5;e6;e7;e8;e9;e10;e11;e12;e13;e14;e15;e16;e17;e18;e19;e20;e21;e22];
F = (10>=x>=0);
F = F + (e>=0);

f = sdpvar(F(find(is(F,'linear'))));
F = F + cut(triu(f*f') >= 0);
sol = optimize(F,obj,sdpsettings('solver','bmibnb','bmibnb.lpreduce',0))

assert(abs(value(obj)--36) <= 1e-6)
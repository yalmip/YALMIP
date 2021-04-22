function tests = test_robust_1
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear')
sdpvar x w;
F = (x + w <= 1) + (-0.5 <= w <= 0.5)
sol = solverobust(F,-x,[],w);

assert(abs(value(x) - 1/2) <=  1e-5)
assert(sol.problem == 0)

F = (abs(x) + w <= 1) + (-0.5 <= w <= 0.5)
sol = solverobust(F,-x,[],w);
assert(abs(value(x) - 1/2) <=  1e-5)
assert(sol.problem == 0)

F = (abs(x) + w <= 1) + (-2 <= w <= 2)
sol = solverobust(F,-x,[],w);
assert((sol.problem == 1) | (sol.problem == 12) | (sol.problem == 15))

sol = optimize([uncertain(w),-2<=w<=4,5-x*w>=0],-x)
assert(abs(value(x) - 1.25) <= 1e-5)

sol = optimize([uncertain(w),-2<=w<=4,5-x*w>=0],x)
assert(abs(value(x) - -2.5) <= 1e-5)

F = (abs(x) + w <= 1) + (norm(w,1) <= 0.3);
sol = solverobust(F,-x,[],w);
assert(abs(value(x) - 0.7) <= 1e-5)
assert(sol.problem ==  0)

F = (abs(x) + w <= 1) + (norm(w,2) <= 0.3);
sol = solverobust(F,-x,[],w);
assert(abs(value(x) - 0.7) <= 1e-5)
assert(sol.problem ==  0)


x = sdpvar(3,1);
w = sdpvar(2,1);
F = (norm(x,1) + sum(w) <= 1) + (norm(w,1) <= 0.3);
sol = solverobust(F,-sum(x),[],w);
assert(abs(value(sum(x)) - 0.7) <= 1e-5)
assert(sol.problem == 0)

x = sdpvar(3,1);
w = sdpvar(2,1);
F = (norm(x,1) + sum(w) <= 1) + (norm(w,1) <= 0.3) + (uncertain(w))
sol = optimize(F,-sum(x));
assert(abs(value(sum(x)) - 0.7) <= 1e-5)
assert(sol.problem == 0)

x = sdpvar(3,1);
w = sdpvar(2,1);
F = (norm(x,1) + sum(w) <= 1) + (w'*w <= 0.3^2);
sol = solverobust(F,-sum(x),[],w);
assert(abs(value(sum(x)) - 0.57574) <= 1e-5)
assert(sol.problem == 0)

x = sdpvar(3,1);
w = sdpvar(2,1);
F = (norm(x,1) + sum(w) <= 1) + (norm(w,2)<=0.3);
sol = solverobust(F,-sum(x),[],w);
assert(abs(value(sum(x)) - 0.57574) <= 1e-5)
assert(sol.problem == 0)

% multiple nonlinear operators in uncertainty descr.
F = (norm(x,1) + sum(w) <= 1) + (norm(w,1) <= 0.3) + (norm(w,inf) <= 0.1);
sol = solverobust(F,-sum(x),[],w);
assert(abs(value(sum(x)) - 0.8) <= 1e-5)
assert(sol.problem == 0)

% Mixed nonlinear operator
sdpvar x w
sol = solverobust((norm(x+w,1)  <= 1) + (norm(w,2) <= 0.2),-sum(x),[],w)
assert(abs(value(sum(x)) - 0.8) <= 1e-5)
assert(sol.problem == 0)

% Mixed nonlinear operator, only part of w in constraints
%sdpvar x w(5,1)
%sol = solverobust((norm(x+w(1),1)  <= 1) + (norm(w,2) <= 0.2),-sum(x),[],w)
%mbg_asserttolequal(value(sum(x)), 0.8, 1e-5);
%assert(sol.problem == 0)
% 
% yalmip('clear')
% 
% x = sdpvar(3,1);
% w = sdpvar(2,1);
% 
% F = (norm(x,1) + norm(w,1) <= 1) + (norm(w,1) <= 0.3) + (norm(w,inf) <= 0.1);
% sol = solverobust(F,-sum(x),[],w);
% mbg_asserttolequal(value(sum(x)), 0.8, 1e-5);
% assert(sol.problem == 0)

% sdpvar x w(2,1)
% F = (x+sum(w) <= 1) + ([1 x;x 2] > 0);
% W = (w'*w <= 1/2);
% objective = (x+1)'*(x+1) + x*norm(w,1);
% sol = solverobust(F + W,objective,[],w);
% mbg_asserttolequal(value(x), -1, 1e-4);

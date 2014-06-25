function test_robust_1

yalmip('clear')
sdpvar x w;
F = set(x + w <= 1) + set(-0.5 <= w <= 0.5)
sol = solverobust(F,-x,[],w);

mbg_asserttolequal(double(x), 1/2, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

F = set(abs(x) + w <= 1) + set(-0.5 <= w <= 0.5)
sol = solverobust(F,-x,[],w);
mbg_asserttolequal(double(x), 1/2, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

F = set(abs(x) + w <= 1) + set(-2 <= w <= 2)
sol = solverobust(F,-x,[],w);
mbg_asserttolequal(double((sol.problem == 1) | (sol.problem == 12) | (sol.problem == 15)), 1, 1e-5);

sol = solvesdp([uncertain(w),-2<=w<=4,5-x*w>=0],-x)
mbg_asserttolequal(double(x), 1.25, 1e-5);

sol = solvesdp([uncertain(w),-2<=w<=4,5-x*w>=0],x)
mbg_asserttolequal(double(x), -2.5, 1e-5);

F = set(abs(x) + w <= 1) + set(norm(w,1) <= 0.3);
sol = solverobust(F,-x,[],w);
mbg_asserttolequal(double(x), 0.7, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

F = set(abs(x) + w <= 1) + set(norm(w,2) <= 0.3);
sol = solverobust(F,-x,[],w);
mbg_asserttolequal(double(x), 0.7, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);


x = sdpvar(3,1);
w = sdpvar(2,1);
F = set(norm(x,1) + sum(w) <= 1) + set(norm(w,1) <= 0.3);
sol = solverobust(F,-sum(x),[],w);
mbg_asserttolequal(double(sum(x)), 0.7, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

x = sdpvar(3,1);
w = sdpvar(2,1);
F = set(norm(x,1) + sum(w) <= 1) + set(norm(w,1) <= 0.3) + set(uncertain(w))
sol = solvesdp(F,-sum(x));
mbg_asserttolequal(double(sum(x)), 0.7, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

x = sdpvar(3,1);
w = sdpvar(2,1);
F = set(norm(x,1) + sum(w) <= 1) + set(w'*w <= 0.3^2);
sol = solverobust(F,-sum(x),[],w);
mbg_asserttolequal(double(sum(x)), 0.57574, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

x = sdpvar(3,1);
w = sdpvar(2,1);
F = set(norm(x,1) + sum(w) <= 1) + set(norm(w,2)<=0.3);
sol = solverobust(F,-sum(x),[],w);
mbg_asserttolequal(double(sum(x)), 0.57574, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

% multiple nonlinear operators in uncertainty descr.
F = set(norm(x,1) + sum(w) <= 1) + set(norm(w,1) <= 0.3) + set(norm(w,inf) <= 0.1);
sol = solverobust(F,-sum(x),[],w);
mbg_asserttolequal(double(sum(x)), 0.8, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

% Mixed nonlinear operator
sdpvar x w
sol = solverobust(set(norm(x+w,1)  <= 1) + set(norm(w,2) <= 0.2),-sum(x),[],w)
mbg_asserttolequal(double(sum(x)), 0.8, 1e-5);
mbg_asserttolequal(sol.problem, 0, 1e-5);

% Mixed nonlinear operator, only part of w in constraints
%sdpvar x w(5,1)
%sol = solverobust(set(norm(x+w(1),1)  <= 1) + set(norm(w,2) <= 0.2),-sum(x),[],w)
%mbg_asserttolequal(double(sum(x)), 0.8, 1e-5);
%mbg_asserttolequal(sol.problem, 0, 1e-5);
% 
% yalmip('clear')
% 
% x = sdpvar(3,1);
% w = sdpvar(2,1);
% 
% F = set(norm(x,1) + norm(w,1) <= 1) + set(norm(w,1) <= 0.3) + set(norm(w,inf) <= 0.1);
% sol = solverobust(F,-sum(x),[],w);
% mbg_asserttolequal(double(sum(x)), 0.8, 1e-5);
% mbg_asserttolequal(sol.problem, 0, 1e-5);

% sdpvar x w(2,1)
% F = set(x+sum(w) <= 1) + set([1 x;x 2] > 0);
% W = set(w'*w <= 1/2);
% objective = (x+1)'*(x+1) + x*norm(w,1);
% sol = solverobust(F + W,objective,[],w);
% mbg_asserttolequal(double(x), -1, 1e-4);

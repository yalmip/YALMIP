function tests = test_robust_13
tests = functiontests(localfunctions);

function test1(testCase)

% Test product between decision variable and auxilliary variable with w
% dependence in objective. This is not possible 
% Really? fmincon seem to get started now
% yalmip('clear')
% sdpvar x w(2,1)
% F = [x+sum(w) <= 1, [1 x;x 2] >= 0];
% W = [w'*w <= 1/2, uncertain(w)];
% objective = (x+1)'*(x+1) + x*norm(w,1);
% sol = optimize(F + W,objective,sdpsettings('verbose',0));
% testCase.assertTrue(sol.problem == -2);

% The norm variable will be an auxilliary variable, not uncertainty
yalmip('clear')
sdpvar x w(2,1)
F = [x + norm(w,1) <= 1];
W = [w'*w <= 1/2, uncertain(w)];
objective = x;
sol = optimize(F + W,-x,sdpsettings('verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(x)--0.41421) <= 1e-3);

sol = optimize(F + W,-x,sdpsettings('robust.auxreduce','projection','verbose',0));
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(x)) <= 1e-3);

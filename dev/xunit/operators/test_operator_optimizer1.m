function tests = test_operator_optimizer1
tests = functiontests(localfunctions);

function test1(dummy)
% Tests a regression bug that made expandmodel flawed. Basically, when
% optmizer generates the model, it constraints the parametric variables to
% pi. However, these constraints can not be used to tighten the model since
% they are completely artificial.

yalmip('clear')
sdpvar x u r
constraints = [abs(x + u - r)>=0.01];
constraints = [constraints, -5<= x <=5, -1 <= u <= 1];
constraints = [constraints, -5<= r <=5];
objective = u;
controller = optimizer(constraints, objective,sdpsettings('verbose',1),[x;r],u);
[u,sol]=controller{[2;3]}

assert(sol == 0);
assert(abs(u--1) <= 1e-4);

function test_optimizer2

% Tests a regression bug that made expandmodel flawed. Basically, when
% optmizer generates the model, it constraints the parametric variables to
% pi. However, these constraints can not be used to tighten the model since
% they are completely artificial.

yalmip('clear')
sdpvar x
sdpvar z
sdpvar y
P = optimizer([-1 <= y <= 1,x <= 10+z*x],x^2,[],[y;z],x)
[~,err] = P{[0;5]};

mbg_asserttrue(err == 0);

yalmip('clear')
sdpvar x
sdpvar y
sdpvar z
P = optimizer([-1 <= y <= 1,x <= 10+z*x],x^2,[],[y;z],x)
[~,err] = P{[0;5]};

mbg_asserttrue(err == 0);

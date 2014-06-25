function test_operator_optimizer3

yalmip('clear')
sdpvar x
sdpvar z
sdpvar y
P = optimizer([-1 <= y <= 1,y >= z*z, x <= x+x*z],x^2,sdpsettings('solver','+sdpt3'),[y;z],x)

% Should be infeasible 
[~,err] = P{[0;2]};
mbg_asserttrue(err == 1 | err == 12);


P = optimizer([-1 <= y <= 1,[y z;z 1]>=0, x <= x+x*z],x^2,sdpsettings('solver','+sdpt3'),[y;z],x)
[~,err] = P{[0;2]};
mbg_asserttrue(err == 1 | err == 12);


function bmibnb_integer

clear sin % To avoid aby problems from code above 
sdpvar x
strange = @(x) sdpfun(x,'@(x) sin(10*x)+abs(sin(x))+x');
obj = strange(x);
sol=solvesdp((integer(x)) + (-pi <= x <= pi),strange(x),sdpsettings('solver','bmibnb'));

mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(x),-2, 1e-4);
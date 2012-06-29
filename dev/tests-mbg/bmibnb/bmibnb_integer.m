function bmibnb_integer

clear sin % To avoid aby problems from code above 
sdpvar x
strange = @(x) sdpfun(x,'@(x) sin(10*x)+abs(sin(x))+x');
obj = strange(x);
sol=solvesdp(set(integer(x)) + set(-pi <= x <= pi),strange(x),sdpsettings('solver','bmibnb'));

mbg_assertfalse(sol.problem);
mbg_asserttolequal(double(x),-2, 1e-4);
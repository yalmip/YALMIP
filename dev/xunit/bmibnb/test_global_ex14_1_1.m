function ex14_1_1

sdpvar x1 x2 x3 objvar

F = ([]);
F = F + (- x3 + objvar == 0);
F = F + (2*sqr(x2) + 4*x1*x2 - 42*x1 + 4*power(x1,3) - x3    <= 14);
F = F + ((-2*sqr(x2)) - 4*x1*x2 + 42*x1 - 4*power(x1,3) - x3 <= -14);
F = F + (2*sqr(x1) + 4*x1*x2 - 26*x2 + 4*power(x2,3) - x3    <= 22);
F = F + ((-2*sqr(x1)) - 4*x1*x2 + 26*x2 - 4*power(x2,3) - x3 <= -22);
F = F + (5 >= x1 >= -5) + (-5 <= x2 <= 5);

sol = solvesdp(F,objvar,sdpsettings('solver','bmibnb'))

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(objvar), 0, 1e-4);

double([x1 x2 x3])

function y = sqr(x)
y = x*x;
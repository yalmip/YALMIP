function ex14_1_8

yalmip('clear')

sdpvar x1 x2 x3 objvar; 

F = ([]);

F = F + ( - x3 + objvar == 0);

F = F + (  (0.0476666666666666 - 0.0649999999999999*x1)*exp(10*x1/(1 + 0.01*x1)) - x1      - x3 <= 0);

F = F + (  x1 - (0.0476666666666666 - 0.0649999999999999*x1)*exp(10*x1/(1 + 0.01*x1))      - x3 <= 0);

F = F + (  (0.143 + (-0.13*x1) - 0.195*x2)*exp(10*x2/(1 + 0.01*x2)) + x1 - 3*x2 - x3  <= 0);

F = F + (  (-(0.143 + (-0.13*x1) - 0.195*x2)*exp(10*x2/(1 + 0.01*x2))) - x1 + 3*x2 - x3 <= 0);

F = F + (-1 <= [x1 x2 ] <= 1);

sol = solvesdp(F,objvar,sdpsettings('solver','bmibnb','bmibnb.upper','fmincon','allownon',1));
mbg_asserttrue(sol.problem==0);
mbg_asserttolequal(double(objvar), 0, 2e-2);
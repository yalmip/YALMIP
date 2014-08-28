function ex14_1_3

yalmip('clear');
sdpvar x1 x2 x3 objvar

F = ([]);
F = F + ( - x3 + objvar == 0); 
F = F + ( 10000*x1*x2 - x3 <= 1); 
F = F + ( - 10000*x1*x2 - x3 <= -1); 
F = F + ( exp(-x1) + exp(-x2) - x3 <= 1.001);
F = F + ( (-exp(-x1)) - exp(-x2) - x3 <= -1.001); 
F = F + (5.49e-6 <= x1 <= 4.553) + (18.21 >= x2 >= 0.0021961);

sol = solvesdp(F,objvar,sdpsettings('solver','bmibnb','bmibnb.absgaptol',1e-8,'bmibnb.relgaptol',1e-8))
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double([objvar ]),0, 1e-3);

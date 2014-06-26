function ex14_1_3

yalmip('clear');
sdpvar x1 x2 x3 objvar

F = set([]);
F = F + set( - x3 + objvar == 0); 
F = F + set( 10000*x1*x2 - x3 <= 1); 
F = F + set( - 10000*x1*x2 - x3 <= -1); 
F = F + set( exp(-x1) + exp(-x2) - x3 <= 1.001);
F = F + set( (-exp(-x1)) - exp(-x2) - x3 <= -1.001); 
F = F + set(5.49e-6 <= x1 <= 4.553) + set(18.21 >= x2 >= 0.0021961);

sol = solvesdp(F,objvar,sdpsettings('solver','bmibnb','bmibnb.absgaptol',1e-8,'bmibnb.relgaptol',1e-8))
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double([objvar ]),0, 1e-3);

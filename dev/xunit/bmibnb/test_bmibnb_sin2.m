function sin2

% Tests mixed monomial evals

sdpvar x
obj = sin(cos(x.^2).^2).^2 + 0.1*(x-2).^2;
sol = solvesdp([-5 <= x <= 5],obj,sdpsettings('allownon',1,'solver','bmibnb','bmibnb.absgaptol',1e-8,'bmibnb.relgaptol',1e-8))

mbg_asserttolequal(sol.problem,0, 1e-4);
mbg_asserttolequal(double(obj),0.0022650, 1e-4);
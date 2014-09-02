function clearsdpvar(X,value)
% CLEARSDPVAR Clear solution

x_lmi_variables = X.lmi_variables;

sol = sdpvar('getSolution');
keep_these = find(~ismember(sol.variables,x_lmi_variables));
sol.optvar = [sol.optvar(keep_these)];
sol.variables = [sol.variables(keep_these)];
sdpvar('setAllSolution',sol);




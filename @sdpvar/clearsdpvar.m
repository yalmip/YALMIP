function clearsdpvar(X,value)
% CLEARSDPVAR Clear solution

% Author Johan Löfberg 
% $Id: clearsdpvar.m,v 1.2 2004-07-01 11:17:09 johanl Exp $  

x_lmi_variables = X.lmi_variables;

sol = sdpvar('getSolution');
keep_these = find(~ismember(sol.variables,x_lmi_variables));
sol.optvar = [sol.optvar(keep_these)];
sol.variables = [sol.variables(keep_these)];
sdpvar('setAllSolution',sol);




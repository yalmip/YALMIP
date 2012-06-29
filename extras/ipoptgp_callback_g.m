function g = ipoptgp_callback_g(x,model)

% Should be made faster by re-using results from dg computation...

% Compute the nonlinear terms in the constraints
[g,geq] = fmincon_congp(x,model);

% Append with linear constraints
g = [g;geq];
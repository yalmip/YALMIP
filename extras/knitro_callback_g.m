function [g,geq,dg,dgeq] = knitro_callback_g(x,model)

global latest_x_g
global latest_dg
global latest_g
global latest_geq
global latest_dgeq

x = x(:);
if exist('latest_x_g') && isequal(x,latest_x_g)
    % Jacobian was computed earlier
   g = latest_g;
   dg = latest_dg;
   geq = latest_geq;
   dgeq = latest_dgeq;   
else
    % Compute the nonlinear terms in the constraints
    [g,geq,dg,dgeq,xevaled] = fmincon_con(x,model);
    latest_x_g = x;
    latest_g = g;
    latest_dg = dg;
    latest_geq = geq;
    latest_dgeq = dgeq;
end
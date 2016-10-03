function g = ipopt_callback_g(x,model)

global latest_x_g
global latest_G
global latest_g
x = x(:);

% Compute the nonlinear terms in the constraints and Jacobians for later   
[g,geq,dg,dgeq] = fmincon_con_liftlayer(x,model);    

% Append with linear constraints
g = [g;geq];
if ~isempty(model.A)
    g = [g;model.A*x - model.b];
end
if ~isempty(model.Aeq)
    g = [g;model.Aeq*x - model.beq];
end

% Append with Jacobians with linear terms
G = [dg';dgeq'];
if ~isempty(model.A)
    G = [G;model.A];
end
if ~isempty(model.Aeq)
    G = [G;model.Aeq];
end

% Save the Jacobian, and information about for which x it was computed
latest_G = sparse(G);
latest_x_g = x;
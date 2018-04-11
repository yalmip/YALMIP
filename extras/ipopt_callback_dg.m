function G = ipopt_callback_dg(x,model)%returnStructOnly,model)

global latest_x_g
global latest_G
global latest_g

x = x(:);
if isequal(x,latest_x_g)
    % Jacobian was computed already in the call for constraints
    G = latest_G;
else
    % Compute the nonlinear terms in the constraints
    [g,geq,dg,dgeq] = fmincon_con_liftlayer(x,model);

    % Append with linear terms
    if isempty(dg)
        G = dgeq';
    else
        G = [dg';dgeq'];
    end
    
    if ~isempty(model.A)
        G = [G;model.A];
    end
    if ~isempty(model.Aeq)     
        G = [G;model.Aeq];
    end
    G = sparse(G);
end
if model.dense
    G = full(G);
end
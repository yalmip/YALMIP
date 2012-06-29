function G = ipopt_callback_dg(x,model)%returnStructOnly,model)

global latest_x_g
global latest_G
global latest_g

% if returnStructOnly
%     G = sparse(ones(length(model.beq)+length(model.b)+length(model.bnonlineq)+length(model.bnonlinineq),length(model.linearindicies)));
%     return
% end

x = x(:);
if isequal(x,latest_x_g)
    % Jacobian was computed already in the call for constraints
    G = latest_G;
else
    % Compute the nonlinear terms in the constraints
    [g,geq,dg,dgeq] = fmincon_con(x,model);

    % Append with linear terms
    G = [dg';dgeq'];
    if ~isempty(model.A)
        G = [G;model.A];
    end
    if ~isempty(model.Aeq)
        G = [G;model.Aeq];
    end
    G = sparse(G);
end
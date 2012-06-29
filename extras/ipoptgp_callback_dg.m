function G = ipopt_callback_dg(x,returnStructOnly,prob)

if returnStructOnly
    G = sparse(ones(max(prob.map)+size(prob.G,1),size(prob.A,2)));
    return
end

% Compute the nonlinear terms in the constraints
[g,geq,dg,dgeq] = fmincon_congp(x,prob);
G = [dg';dgeq'];
G = sparse(G);
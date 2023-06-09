function G = ipoptgp_callback_dg(x,prob)

% Compute the nonlinear terms in the constraints
[g,geq,dg,dgeq] = fmincon_congp(x,prob);
G = [dg';dgeq'];
G = sparse(G);

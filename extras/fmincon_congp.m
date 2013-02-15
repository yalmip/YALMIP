function [g,geq,dg,dgeq] = fmincon_congp(x,prob)

z = prob.A*x;
expz = exp(z);

if size(prob.B,1)>0
    g = prob.B*expz;
    g = full(log(g));    
    % Should be correct, but it fails for some problems (test_gp_5)
    dg = [];
    %z = prob.A*x;
    %expz = exp(z);
    expz(isinf(expz)) = 1e5;
    T=diag(sparse(1./(prob.B*expz)));
    U=diag(sparse(expz));
    dg = (T*prob.B*(U*prob.A))';
else
    g = [];
    dg = [];
end

if  length(prob.h) > 0
    geq = log(prob.h) + prob.G*x;
    dgeq = prob.G';
else
    geq = [];
    dgeq = [];
end
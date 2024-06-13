function p = complementmodel(p)

% Extend model to include complemented binary
% Used when doing clique analysis
A = p.F_struc(:,2:end);
B = A(:,p.binary_variables);
p.F_struc(:,1+p.binary_variables) = min(0,B);
B = -max(B,0);
p.F_struc = [p.F_struc B];
p.F_struc(:,1) = p.F_struc(:,1) - sum(B,2);
p.binary_variables = [p.binary_variables length(p.c) + p.binary_variables];
p.isbinary = [p.isbinary p.isbinary ];
lb = p.lb;
ub = p.ub;
p.lb = [p.lb;1-ub];
p.ub = [p.ub;1-lb];
p.cliques = [p.cliques;p.cliques*0];
p.c = [p.c;p.c];


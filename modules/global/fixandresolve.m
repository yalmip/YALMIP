function [upperhere,xtempwork] = monotoneSDPHeuristics(p,upper,x,aux1,aux2)
xtempwork = [];
upperhere = inf;

x_fix = x;
xtempwork = [];
upperhere = inf;
% Heuristically set all cardinality groups
for i = 1:length(p.cardinalitygroups)
    s = p.cardinalitygroups{i};
    n = p.cardinalitysize{i};
    z = x_fix(s);
    [~,loc] = sort(z);
    z(loc(end-n+1:end)) = 1;
    z(loc(1:end-n)) = 0;
    x_fix(s) = z;
end
% Not much more we can do than round the rest
x_fix(p.binary_variables) = round(x_fix(p.binary_variables));
x_fix(p.integer_variables) = round(x_fix(p.integer_variables));

% Force all integer variables
p.lb(p.binary_variables) = x_fix(p.binary_variables);
p.ub(p.binary_variables) = x_fix(p.binary_variables);
p.lb(p.integer_variables) = x_fix(p.integer_variables);
p.ub(p.integer_variables) = x_fix(p.integer_variables);

relaxed_p = p;
relaxed_p.integer_variables = [];
relaxed_p.binary_variables = [];

output = bnb_solvelower(p.solver.lower.call,relaxed_p,inf,-inf,[],0,[]);
xtempwork = output.Primal;
xtempwork = setnonlinearvariables(p,xtempwork);
if checkfeasiblefast(p,xtempwork,p.options.bnb.feastol)
    upperhere = computecost(p.f,p.corig,p.Q,xtempwork,p);
end
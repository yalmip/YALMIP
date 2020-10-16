function [upperhere,xtempwork] = monotoneSDPHeuristics(p,upper,x,aux1,aux2)
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

% Try to find a step-length which rednders us feasible
trials = 0;
s = find(p.c);
perturb = zeros(length(p.c),1);
perturb(s) = abs(x_fix(s)) + 1e-4;
upperhere = inf;
steplength = 0;
while trials < 20
    steplength_lower = steplength;
    steplength = .001*2^(trials+1);
    xtemp = x_fix + steplength*perturb;
    xtemp = setnonlinearvariables(p,xtemp);
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)
        upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
        steplength_upper = steplength;
        break
    end
    trials = trials + 1;
end

if trials < 10
    % By increasing step, we found a feasible solution
    % steplength_lower: Infeasible
    % steplength_upper: Feasible
    xtempwork = xtemp;
    for j = 1:10
        steplength = (steplength_lower+steplength_upper)/2;
        xtemp = x_fix + steplength*perturb;
        xtemp = setnonlinearvariables(p,xtemp);
        if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)
            upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
            xtempwork = xtemp;
            steplength_upper = steplength;
        else
            steplength_lower = steplength;
        end
    end
end


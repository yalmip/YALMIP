function [upper,x_min] = sosfixer(p,relaxedsolution,prelaxed)

% Extremely simple heuristic for finding integer
% solutions in problems where all variables are sos1 constrained

% This was the relaxed solution
x = relaxedsolution.Primal;

% find the largest relaxed value in each sos1 group, and fix that value to
% 1
for i = 1:length(prelaxed.sosgroups)    
    [~,loc] = max(relaxedsolution.Primal(prelaxed.sosgroups{i}));
    loc = prelaxed.sosgroups{i}(loc);
    candidates = find(abs(relaxedsolution.Primal(loc)-relaxedsolution.Primal(prelaxed.sosgroups{i}))<1e-5);
    %candidates(ceil(length(candidates)*rand(1)))
    loc = prelaxed.sosgroups{i}(candidates(ceil(length(candidates)*rand(1))));
    if p.lb(loc)==0
        prelaxed.lb(loc)=1;
        prelaxed.ub(setdiff(prelaxed.sosgroups{i},loc))=0;
%         prelaxed.K.f = prelaxed.K.f + 1;        
%         e = zeros(1,size(p.F_struc,2));
%         e(1) = 1;
%         e(1,1+loc) = -1;
%         prelaxed.F_struc = [e;prelaxed.F_struc];
    end
end

prelaxed.ub(p.binary_variables) = min(1,prelaxed.ub(p.binary_variables));
prelaxed.lb(p.binary_variables) = min(0,prelaxed.lb(p.binary_variables));

upper = inf;
x_min = [];
output = feval(p.solver.lower.call,prelaxed);
xtemp1 = output.Primal;
upper1 = computecost(p.f,p.corig,p.Q,xtemp1,p);
if checkfeasiblefast(p,xtemp1,p.options.bnb.feastol)
    x_min = xtemp1;
    upper = upper1;
end
end


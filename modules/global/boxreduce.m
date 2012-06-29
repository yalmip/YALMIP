function [p,feasible] = boxreduce(p,upper,lower,lpsolver,options,xmin);

if options.bmibnb.lpreduce

    vol_start    = prod(p.ub(p.branch_variables)-p.lb(p.branch_variables));
    diag_before  =  sum(p.ub(p.branch_variables)-p.lb(p.branch_variables));

    [pcut,feasible,lower] = lpbmitighten(p,lower,upper,lpsolver,xmin);
    diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));
    iterations = 0;
    while (diag_after/(1e-18+diag_before) < .75) & feasible & iterations < 4
        %improvethese = ~[(p.lb==pcut.lb)&(p.ub==pcut.ub)];
        [pcut,feasible,lower] = lpbmitighten(pcut,lower,upper,lpsolver,xmin);
        diag_before = diag_after;
        diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));
        iterations = iterations + 1;
    end

    % Clean up...
    for i = 1:length(pcut.lb)
        if (pcut.lb(i)>pcut.ub(i)) & (pcut.lb-pcut.ub < 1e-3)
            pcut.lb(i)=pcut.ub(i);
            pcut = updatemonomialbounds(pcut);
        end
    end
    p.lb = pcut.lb;
    p.ub = pcut.ub;

    p.lb(p.lb<-1e12) = -inf;
    p.ub(p.ub>1e12) = inf;
    p.evalMap = pcut.evalMap;
else   
    feasible = 1;
end

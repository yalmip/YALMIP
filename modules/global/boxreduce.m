function [p,feasible] = boxreduce(p,upper,lower,lpsolver,options,xmin);

if options.bmibnb.lpreduce
    
    improvethese = zeros(length(p.c),1);
    improvethese(p.linears) = 1;
    
    vol_start    = prod(p.ub(p.branch_variables)-p.lb(p.branch_variables));
    diag_before  =  sum(p.ub(p.branch_variables)-p.lb(p.branch_variables));
    span_before = p.ub - p.lb;
    
    [pcut,feasible,lower] = lpbmitighten(p,lower,upper,lpsolver,xmin,improvethese);
    p.counter = pcut.counter;
    diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));
    span_after = pcut.ub - pcut.lb;
    
    % Variables with gap smaller than desired precision is no longer bound
    % propagated
    improvethese(find(abs(span_after) <=  p.options.bmibnb.vartol)) = 0;
    
    not_decreasing = find(span_after >= 0.8*span_before);
    improvethese(not_decreasing) = max(0,improvethese(not_decreasing)-1);
    
    
    iterations = 0;    
%    while (diag_after/(1e-18+diag_before) < .75) & feasible & iterations < 4
    while any(improvethese) & feasible & iterations < 8                
        span_before = span_after;        
        [pcut,feasible,lower] = lpbmitighten(pcut,lower,upper,lpsolver,xmin,improvethese);
        p.counter = pcut.counter;
        span_after = pcut.ub - pcut.lb;
        improvethese(find(abs(span_after) <=  p.options.bmibnb.vartol)) = 0;    
        not_decreasing = find(span_after >= 0.95*span_before);
        improvethese(not_decreasing) = max(0,improvethese(not_decreasing)-1);
        %diag_before = diag_after;
        %diag_after = sum(pcut.ub(p.branch_variables)-pcut.lb(p.branch_variables));
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

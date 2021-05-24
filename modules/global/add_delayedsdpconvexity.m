function p = add_delayedsdpconvexity(p,spliton)

if nargin == 1
    for i = p.linears
        p = checkSingle(p,i);    
    end
else
    p = checkSingle(p,spliton);
end

function p = checkSingle(p,spliton)
if any(p.hiddendelayedconvex.variable)
    loc = find(spliton == p.hiddendelayedconvex.variable(1,:));
    if any(loc)
        % x >= u or x <= l
        u = p.hiddendelayedconvex.variable(2,loc);
        l = p.hiddendelayedconvex.variable(3,loc);
        if p.ub(spliton) <= l
            p.delayedconvex = mergeNumericalModels(p.delayedconvex,p.hiddendelayedconvex.modeldown{loc});
            p.hiddendelayedconvex.variable(3,loc)=nan;
        end
        if p.lb(spliton) >= u
            p.delayedconvex = mergeNumericalModels(p.delayedconvex,p.hiddendelayedconvex.modelup{loc});
            p.hiddendelayedconvex.variable(2,loc)=nan;
        end
    end
    a=isnan(p.hiddendelayedconvex.variable(2,loc));
    b=isnan(p.hiddendelayedconvex.variable(3,loc));
    if any(a) && any(b)
        p.hiddendelayedconvex.variable(1,loc)=0;
    end
end
function [upper,x_min] = upper_from_sdpextension(p,xtemp,upper,x_min)
if p.sdpextendable && p.lb(p.noninteger_variables)<p.ub(p.noninteger_variables)
    [xtemp,fail] = sdpextendsolution(p,xtemp);
    if ~fail
        if ~isnan(xtemp(p.noninteger_variables))
            xtemp(p.noninteger_variables) = max(xtemp(p.noninteger_variables),p.lb(p.noninteger_variables));
            xtemp(p.noninteger_variables) = min(xtemp(p.noninteger_variables),p.ub(p.noninteger_variables));
            upperhere = computecost(p.f,p.c,p.Q,xtemp,p);
            if upperhere < upper && checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
                x_min = xtemp;
                upper = upperhere;
            end
        end
    end
end
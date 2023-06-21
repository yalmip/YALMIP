function [local_upper,x_min,covers] = rounder(prelaxed,upper,x,p,relaxedoutput,lower,covers,dummy)
% Extremely simple heuristic for finding integer solutions.

% This was the relaxed solution
x = relaxedoutput.Primal;

% Assume we fail
local_upper = inf;
x_min = x;

% These should be integer
intvars = p.integral_variables(:);
convars = p.noninteger_variables;

% Clean up things that we can consider integer
close = find(abs(x(intvars)-round(x(intvars)))<=p.options.bnb.inttol);
x(intvars(close)) = round(x(intvars(close)));

xtemp = x;
xtemp(intvars) = round(x(intvars));

if ismember('shifted round',p.options.bnb.rounding)

    % Round, update nonlinear terms, and compute feasibility
    oldxtemp = [];
    for tt = -.4:0.1:0.4
        xtemp = x; 
        xtemp(intvars) = xtemp(intvars)+tt;       
        xtest = xtemp;
        xtemp = min(xtemp,p.ub);
        xtemp = max(xtemp,p.lb);
        xtemp(intvars) = round(xtemp(intvars));
        xtemp = min(xtemp,p.ub);
        xtemp = max(xtemp,p.lb);               
        xtemp = fix_binary_products(p,xtemp);        
        xtemp = fix_downforce(p,xtemp);
        xtemp = fix_cardinality(p,xtemp,xtest);
        xtemp = fix_semivar(p,xtemp);        
        xtemp = setnonlinearvariables(p,xtemp);
              
        if ~isequal(xtemp,oldxtemp)
            oldxtemp=xtemp;
            upperhere = computecost(p.f,p.c,p.Q,xtemp,p);
            if upperhere < upper && upperhere >= lower
                if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)
                    x_min = xtemp;
                    local_upper = upperhere;
                    upper=local_upper;
                   
                elseif ~isequal(xtemp(intvars),x_min(intvars))
                    % Check for common SDP case such as maximizing smallest eigenvalue
                    % or minimizing largest.
                    % With x fixed, smallest t can be computed by gevp                    
                   % [xtemp,fail] = sdpextendsolution(p,xtemp);
                   [upperhere,xhere] = upper_from_sdpextension(p,xtemp,upper,x_min);
                   if upperhere<upper
                       x_min = xhere;
                       local_upper = upperhere;
                       upper = local_upper;
                   end
%                     if ~fail
%                         if ~isnan(xtemp(convars))
%                             xtemp(convars) = max(xtemp(convars),p.lb(convars));
%                             xtemp(convars) = min(xtemp(convars),p.ub(convars));
%                             upperhere = computecost(p.f,p.c,p.Q,xtemp,p);
%                             if upperhere < upper && checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
%                                 x_min = xtemp;
%                                 local_upper = upperhere;
%                                 upper = local_upper;
%                             end
%                         end
%                     end                    
                end
            end
        end
    end
end

function x = fix_semivar(p,x)
for i = 1:length(p.semicont_variables)
    j = p.semicont_variables(i);
    if x(j)>= p.semibounds.lb(i) && x(j)<= p.semibounds.ub(i)
        % OK
    elseif x(j)==0
        % OK
    else
        s = [abs(x(j)-0); abs(x(j)-p.semibounds.lb(i));abs(x(j)-p.semibounds.ub(i))];
        [dummy,index] = min(s);
        switch index
            case 1
                x(j) = 0;
            case 2
                x(j) = p.semibounds.lb(i);
            case 3
                x(j) = p.semibounds.lb(i);
        end
    end
end
function [upper,x_min] = rounder(p,relaxedsolution,prelaxed)

% Extremely simple heuristic for finding integer
% solutions.
%
% Rounds up and down, fixes etc.

% This was the relaxed solution
x = relaxedsolution.Primal;

% Assume we fail
upper = inf;
x_min = x;

% These should be integer
intvars = [p.integer_variables(:);p.binary_variables(:)];

if ismember('shifted ceil',p.options.bnb.rounding)
    % Round, update nonlinear terms, and compute feasibility
    for tt = logspace(0,-4,4)
        
        f = x(intvars)-floor(x(intvars));
        xtemp = x;xtemp(intvars) = round(xtemp(intvars));
        xtemp(intvars(f > tt)) = ceil(x(intvars(f > tt)));
        xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
        xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
        xtemp = fix_semivar(p,xtemp);
        xtemp = setnonlinearvariables(p,xtemp);
        if isfield(p.options,'plottruss')
            if p.options.plottruss
                plottruss(4,'Shifted ceil',p,xtemp);
            end
        end
        upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
        if upperhere < upper &  checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
            x_min = xtemp;
            upper =upperhere;            
        end
    end
end

if ismember('shifted round',p.options.bnb.rounding)
    % Round, update nonlinear terms, and compute feasibility
    for tt = 0:0.05:0.45
        xtemp = x;xtemp(intvars) = round(xtemp(intvars)+tt);
        xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
        xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
        xtemp = fix_semivar(p,xtemp);
        xtemp = setnonlinearvariables(p,xtemp);
        if isfield(p.options,'plottruss')
            if p.options.plottruss
                plottruss(2,'Shifted round',p,xtemp);
            end
        end
        upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
        if upperhere < upper &  checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
            x_min = xtemp;
            upper =upperhere;%p.f+x_min'*p.Q*x_min + p.corig'*x_min;
            return      
        end
    end
end

if length(prelaxed.sosgroups)>0
    xtemp = x;
    for i = 1:length(prelaxed.sosgroups)
        xi = x(prelaxed.sosgroups{1});
        [j,loc] = max(xi);
        xtemp(prelaxed.sosgroups{i}) = 0;
        xtemp(prelaxed.sosgroups{i}(loc)) = 1;
    end
    xtemp = setnonlinearvariables(p,xtemp);
    upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
    if upperhere < upper &  checkfeasiblefast(p,xtemp,p.options.bnb.feastol)
        x_min = xtemp;
        upper =upperhere;
        return
    end
end

if upper<inf
    return
end


if ismember('round',p.options.bnb.rounding)
    
    % Round, update nonlinear terms, and compute feasibility
    xtemp = x;xtemp(intvars) = round(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp);
    if isfield(p.options,'plottruss')
        if p.options.plottruss
            subplot(2,2,2);
            cla
            title('Rounded node')
            pic(p.options.truss,xtemp);
            drawnow
        end
    end
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = computecost(p.f,p.corig,p.Q,x_min,p);%p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return    
    end
end

if ismember('fix',p.options.bnb.rounding)
    % Do same using fix instead
    xtemp = x;xtemp(intvars) = fix(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp);
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = computecost(p.f,p.corig,p.Q,x_min,p);%upper = p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return
    end
end

if ismember('ceil',p.options.bnb.rounding)
    % ...or ceil
    xtemp = x;xtemp(intvars) = ceil(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp);
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = computecost(p.f,p.corig,p.Q,x_min,p);%upper = p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return   
     end
end

if ismember('floor',p.options.bnb.rounding)
    % or floor
    xtemp = x;xtemp(intvars) = floor(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp);
    if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = computecost(p.f,p.corig,p.Q,x_min,p);%upper = p.f+x_min'*p.Q*x_min + p.corig'*x_min;
        return    
     end
end

function x = fix_semivar(p,x);
for i = 1:length(p.semicont_variables)
    j = p.semicont_variables(i);
    if x(j)>= p.semibounds.lb(i) & x(j)<= p.semibounds.ub(i)
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
function [upper,x_min] = rounder(p,relaxedsolution,prelaxed)

% Extremely simple heuristic for finding integer solutions.

% This was the relaxed solution
x = relaxedsolution.Primal;

% Assume we fail
upper = inf;
x_min = x;

% These should be integer
intvars = [p.integer_variables(:);p.binary_variables(:)];
convars = p.noninteger_variables;

if ismember('shifted round',p.options.bnb.rounding)
    % Pre-extract...    
    if length(convars)==1 && length(p.K.s)==1 && p.K.s(1) > 0
        H = p.F_struc(1+p.K.l+p.K.f:end,:);
        H0 = reshape(H(:,1),p.K.s(1),p.K.s(1));if nnz(H0)/numel(H0)>0.5;H0 = full(H0);end
        Hx = reshape(H(:,1+convars),p.K.s(1),p.K.s(1));if nnz(Hx)/numel(Hx)>0.5;Hx = full(Hx);end
        Hz = H(:,1 + intvars);if nnz(Hz)/numel(Hz)>0.5;Hz = full(Hz);end
    end
    % Round, update nonlinear terms, and compute feasibility
    for tt = -.5:0.1:0.5
        xtemp = x;xtemp(intvars) = round(xtemp(intvars)+tt);
        xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
        xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
        xtemp = fix_semivar(p,xtemp);  
        xtemp = fix_atmost(p,xtemp,x);
        xtemp = setnonlinearvariables(p,xtemp);        
        upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
        if upperhere < upper
            if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
                x_min = xtemp;
                upper =upperhere;               
            else
                % Check for common SDP case such as maximizing smallest eigenvalue
                % or minimizing largest.                  
                % With x fixed, smallest t can be computed by gevp
                % TODO: Support and loop over several LMIs
                if length(convars) == 1 && length(p.K.s)==1 && p.K.s(1)>0                    
                    Hy = H0 + reshape(Hz*xtemp(intvars),p.K.s(1),p.K.s(1));
                    s = eig(full(Hx),full(Hy));
                    s(isinf(s))=[];
                    s(isnan(s))=[];
                    if any(s)                             
                        xtemp(convars) = min(-1./s(s~=0));                        
                        if ~isnan(xtemp(convars))
                            xtemp(convars) = max(xtemp(convars),p.lb(convars));
                            xtemp(convars) = min(xtemp(convars),p.ub(convars));
                            upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
                            if upperhere < upper && checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
                                x_min = xtemp;
                                upper = upperhere;                          
                            end
                        end                    
                    end
                end
            end
        end
    end
end

if length(prelaxed.sosgroups)>0
    try
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
        if all(xtemp(p.binary_variables) == fix(xtemp(p.binary_variables))) && all(xtemp(p.integer_variables) == fix(xtemp(p.integer_variables)))
            x_min = xtemp;
            upper =upperhere;
        end
    end
    catch
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

function xtemp = fix_atmost(p,xtemp,x);

would = zeros(1,length(x));
for i = 1:length(p.atmost.groups)   
    k = p.atmost.groups{i};    
    if nnz(xtemp(k))> p.atmost.bounds(i);
        n_should_be_zero = length(p.atmost.groups{i}) - p.atmost.bounds(i);
        [y,loc] = sort(abs(x(k)));
        xtemp(k(loc(1:n_should_be_zero))) = 0;
    end
end
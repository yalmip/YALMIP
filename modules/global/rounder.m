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
        upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
        if upperhere < upper 
            if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
                x_min = xtemp;
                upper =upperhere;                                       
            end
        end
    end
end

if ismember('shifted round',p.options.bnb.rounding)
    % Pre-extract...
    if nnz(p.Q)==0 && nnz(p.corig)==1 && length(p.K.s)==1 && p.K.s(1)>0
        k = find(p.corig);
        if ~ismember(k,intvars)
            other = setdiff(1:length(p.corig),k);
            H = p.F_struc(1+p.K.l+p.K.f:end,:);
            H0 = reshape(H(:,1),p.K.s(1),p.K.s(1));
            Hx = reshape(H(:,1+k),p.K.s(1),p.K.s(1));
            Hz = H(:,1 + other);
        end
    end
    % Round, update nonlinear terms, and compute feasibility
    for tt = -.5:0.1:0.5
        xtemp = x;xtemp(intvars) = round(xtemp(intvars)+tt);
        xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
        xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
        xtemp = fix_semivar(p,xtemp);
      %  xtemp = fix_atmost(p,xtemp,x);
        xtemp = setnonlinearvariables(p,xtemp);  
        if nnz(xtemp(intvars)) > p.cardinality.upper
            [sorted,loc] = sort(abs(x(intvars)));   
            xtemp(intvars(loc(1:(length(intvars)-p.cardinality.upper))))=0;
        end
        upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
        if upperhere < upper
            if checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
                x_min = xtemp;
                upper =upperhere;%p.f+x_min'*p.Q*x_min + p.corig'*x_min;
               % return
            else
                % Check for common SDP case such as maximizing smallest eigenvalue
                % or minimizing largest.
                % We are lookng for 1 variable objective, and that variable
                % only enters a single matrix constraint                
                % With x fixed, smallest t can be computed by gevp
                if nnz(p.Q)==0 && nnz(p.corig)==1 && length(p.K.s)==1 && p.K.s(1)>0
                    k = find(p.corig);
                    if ~ismember(k,intvars)                        
                      %  other = setdiff(1:length(p.corig),k);
                      %  H = p.F_struc(1+p.K.l+p.K.f:end,:);
                      %  H0 = reshape(H(:,1),p.K.s(1),p.K.s(1));
                      %  Hx = reshape(H(:,1+k),p.K.s(1),p.K.s(1));
                      %  Hy = H0 + reshape(H(:,1 + other)*xtemp(other),p.K.s(1),p.K.s(1));
                        Hy = H0 + reshape(Hz*xtemp(other),p.K.s(1),p.K.s(1));
                        xtemp(k) = min(-1./eig(full(Hx),full(Hy)));
                        upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
                        if upperhere < upper && checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
                            x_min = xtemp;
                            upper =upperhere;
                        %    return
                        end
                    end
                end
            end
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
        if all(xtemp(p.binary_variables) == fix(xtemp(p.binary_variables))) && all(xtemp(p.integer_variables) == fix(xtemp(p.integer_variables)))
            x_min = xtemp;
            upper =upperhere;
            return
            end
    end
end

if ismember('round',p.options.bnb.rounding)    
    % Round, update nonlinear terms, and compute feasibility
    xtemp = x;xtemp(intvars) = round(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp); 
    upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
    if upperhere < upper && checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = upperhere;
    end
end

if ismember('fix',p.options.bnb.rounding)
    % Do same using fix instead
    xtemp = x;xtemp(intvars) = fix(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp);
    upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
    if upperhere < upper & checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = upperhere;
    end
end

if ismember('ceil',p.options.bnb.rounding)
    % ...or ceil
    xtemp = x;xtemp(intvars) = ceil(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp);
    upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
    if upperhere < upper && checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = upperhere;
    end
end

if ismember('floor',p.options.bnb.rounding)
    % or floor
    xtemp = x;xtemp(intvars) = floor(xtemp(intvars));
    xtemp(p.binary_variables(:)) = min(1,xtemp(p.binary_variables(:)));
    xtemp(p.binary_variables(:)) = max(0,xtemp(p.binary_variables(:)));
    xtemp = fix_semivar(p,xtemp);
    xtemp = setnonlinearvariables(p,xtemp);
    upperhere = computecost(p.f,p.corig,p.Q,xtemp,p);
    if upperhere < upper && checkfeasiblefast(p,xtemp,p.options.bnb.feastol)%if res>-p.options.bnb.feastol
        x_min = xtemp;
        upper = upperhere;
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
        would(k(loc(1:n_should_be_zero))) =  would(k(loc(1:n_should_be_zero))) + 1;        
    end
end

for i = 1:length(p.atmost.groups)
    k = p.atmost.groups{i};
    if nnz(xtemp(k))> p.atmost.bounds(i);
        adjusted = 1;
        n_should_be_zero = length(p.atmost.groups{i}) - p.atmost.bounds(i);
        [y,loc] = sort(-would(k)+abs(x(k))');
        xtemp(k(loc(1:n_should_be_zero))) = 0;
        x(k(loc(1:n_should_be_zero))) = 0;
        would(k(loc(n_should_be_zero+1:end))) =  would(k(loc(n_should_be_zero+1:end)))-1;           
    end
end

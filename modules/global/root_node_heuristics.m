function [x_min,upper] = root_node_heuristics(p,x_min,upper)

% Structure when when we can extend a single continuous
intvars = [p.integer_variables(:);p.binary_variables(:)];
convars = p.noninteger_variables;
if length(convars)==1 && length(p.K.s)==1 && any(p.K.s)
    sdp.H = p.F_struc(1+p.K.l+p.K.f:end,:);
    sdp.H0 = reshape(sdp.H(:,1),p.K.s(1),p.K.s(1));if nnz(sdp.H0)/numel(sdp.H0)>0.5;sdp.H0 = full(sdp.H0);end
    sdp.Hx = reshape(sdp.H(:,1+convars),p.K.s(1),p.K.s(1));if nnz(sdp.Hx)/numel(sdp.Hx)>0.5;sdp.Hx = full(sdp.Hx);end
    sdp.Hz = sdp.H(:,1 + intvars);if nnz(sdp.Hz)/numel(sdp.Hz)>0.5;sdp.Hz = full(sdp.Hz);end
else
    sdp = [];
end

% Guess 0, except those fixed to something else
x = zeros(length(p.c),1);
x(p.lb == p.ub) = p.lb(p.lb==p.ub);
x = set_binary_products(x,p);
upperhere = computecost(p.f,p.corig,p.Q,x,p);
if upperhere < upper
    if checkfeasiblefast(p,x,p.options.bnb.feastol)
        x_min = x;
        upper = upperhere;
    elseif ~isempty(sdp)
        [upper,x_min] = extendSDP(p,x,sdp,convars,intvars,upper,x_min);
    end
end

% Guess 1, except those fixed to something else
x = zeros(length(p.c),1);
x(p.binary_variables) = 1;
x(p.lb == p.ub) = p.lb(p.lb==p.ub);
x = set_binary_products(x,p);
upperhere = computecost(p.f,p.corig,p.Q,x,p);
if upperhere < upper
    if checkfeasiblefast(p,x,p.options.bnb.feastol)
        x_min = x;
        upper = upper_here;
    elseif ~isempty(sdp)
        [upper,x_min] = extendSDP(p,x,sdp,convars,intvars,upper,x_min);
    end
end

% Guess 0, except 1 n every cardinality knapsack where we pick cheapest
% Super greedy, we start by sorting knapsacks too
if ~isempty(p.knapsack.variables)
    x = zeros(length(p.c),1);
    x(p.lb == p.ub) = p.lb(p.lb==p.ub);
    for i = 1:length(p.knapsack.variables)
        c(i) = min(p.c(p.knapsack.variables{i}));
    end
    [~,cloc] = sort(c,'ascend');
    for i = cloc
        if p.knapsack.type(i)
            if nnz(x(p.knapsack.variables{i})) == 0
                [~,loc] = sort(p.c(p.knapsack.variables{i}),'ascend');
                x(p.knapsack.variables{i}(loc(1))) = 1;
            end
        end
        x = set_binary_products(x,p);
        upperhere = computecost(p.f,p.corig,p.Q,x,p);
        if upperhere < upper && ~isempty(sdp)
            [upper,x_min] = extendSDP(p,x,sdp,convars,intvars,upper,x_min);
        end
    end
    upperhere = computecost(p.f,p.corig,p.Q,x,p);
    if upperhere < upper && isempty(sdp) && checkfeasiblefast(p,x,p.options.bnb.feastol)
        x_min = x;
        upper = upperhere;
    elseif upperhere < upper && ~isempty(sdp)
        [upper,x_min] = extendSDP(p,x,sdp,convars,intvars,upper,x_min);
    end
    % Guess 0, except 1 n every cardinality knapsack where we pick most
    % expensive
      x = zeros(length(p.c),1);
    x(p.lb == p.ub) = p.lb(p.lb==p.ub);
    [~,cloc] = sort(c,'descend');
    for i = cloc
        if p.knapsack.type(i)
            if nnz(x(p.knapsack.variables{i})) == 0
                [~,loc] = sort(p.c(p.knapsack.variables{i}),'descend');
                x(p.knapsack.variables{i}(loc(1))) = 1;
            end
        end
        x = set_binary_products(x,p);
        upperhere = computecost(p.f,p.corig,p.Q,x,p);
        if upperhere < upper && ~isempty(sdp)
            [upper,x_min] = extendSDP(p,x,sdp,convars,intvars,upper,x_min);
        end
    end
    upperhere = computecost(p.f,p.corig,p.Q,x,p);
    if upperhere < upper && isempty(sdp) && checkfeasiblefast(p,x,p.options.bnb.feastol)
        x_min = x;
        upper = upperhere;
    elseif upperhere < upper && ~isempty(sdp)
        [upper,x_min] = extendSDP(p,x,sdp,convars,intvars,upper,x_min);
    end
end

function [upper,x_min] = extendSDP(p,xtemp,sdp,convars,intvars,upper,x_min)

Hy = sdp.H0 + reshape(sdp.Hz*xtemp(intvars),p.K.s(1),p.K.s(1));
s = eig(full(sdp.Hx),full(Hy));
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

function x = set_binary_products(x,p)
if ~isempty(p.binaryProduct)
    x(p.binaryProduct(:,1)) = prod(x(p.binaryProduct(:,2:3)),2);
end
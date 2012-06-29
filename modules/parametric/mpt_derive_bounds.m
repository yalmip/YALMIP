function [Matrices,infeasible]  = derive_bounds(Matrices,options)

A = [Matrices.G -Matrices.E];
b = [Matrices.W];
Aeq = [Matrices.Aeq Matrices.Beq];
beq = [Matrices.beq];
lb = Matrices.lb;
ub = Matrices.ub;

oldlb = lb;
oldub = ub;

fixed = find(lb == ub);
notfixed = setdiff(1:length(ub),fixed);

if ~isempty(fixed)
    AA = A;
    bb = b;
    nn =     lb(fixed);
    b = b-A(:,fixed)*lb(fixed);
    A(:,fixed) = [];

    lb = lb(notfixed);
    ub = ub(notfixed);
end

if ~isempty(fixed) & ~isempty(Aeq)
    beq = beq-Aeq(:,fixed)*nn;
    Aeq(:,fixed) = [];
end

n = size(A,2);
In = eye(n);

% Detect slacks and avoid maximizing them (they are unbounded)
trivially_unbounded_above = zeros(n,1);
trivially_unbounded_above(find(sum(A,1) == -sum(abs(A),1) & isinf(ub)')) = 1;
trivially_unbounded_below = zeros(n,1);
trivially_unbounded_below(find(sum(A,1) == sum(abs(A),1) & isinf(lb)')) = 1;
if ~isempty(Aeq)
    trivially_unbounded_below(find(any(Matrices.Aeq,1))) = 0;
    trivially_unbounded_above(find(any(Matrices.Aeq,1))) = 0;
end

reset_to_minus_inf = find(isinf(ub(find(trivially_unbounded_below))));
reset_to_inf = find(isinf(ub(find(trivially_unbounded_above))));

lb=max(-1e7,lb);
ub=min(1e7,ub);
u = ub;
l = lb;

for i=1:n,

    if ~trivially_unbounded_above(i)
        f=-In(:,i);

        finite_lb = find(~isinf(l));
        finite_ub = find(~isinf(u));
        Abounds = [A;In(finite_ub,:);-In(finite_lb,:)];
        bbounds = [b;u(finite_ub,:);-l(finite_lb,:)];

        Abounds = [A;eye(n);-eye(n)];
        bbounds = [b;u;-l];

        [x,fval,lambda,exitflag,how]=mpt_solveLP(f(:)',Abounds,bbounds,[],[],[],options.mpt.lpsolver);
        if isequal(how,'infeasible')
            Matrices.lb = repmat(inf,length(oldlb),1);
            Matrices.ub = repmat(-inf,length(oldlb),1);
            infeasible = 1;
            return
        end
        u(i,1)=min(u(i,1),x(i));
    end

    if ~trivially_unbounded_below(i)
        f=In(:,i);
        finite_lb = find(~isinf(l));
        finite_ub = find(~isinf(u));
        Abounds = [A;In(finite_ub,:);-In(finite_lb,:)];
        bbounds = [b;u(finite_ub,:);-l(finite_lb,:)];

        Abounds = [A;eye(n);-eye(n)];
        bbounds = [b;u;-l];

        [x,fval,lambda,exitflag,how]=mpt_solveLP(f(:)',Abounds,bbounds,[],[],[],options.mpt.lpsolver);
        if isequal(how,'infeasible')
            Matrices.lb = repmat(inf,length(oldlb),1);
            Matrices.ub = repmat(-inf,length(oldlb),1);
            infeasible = 1;
            return
        end
        l(i,1)=max(lb(i),x(i));
    end

end

ll = oldlb;
if ~isempty(notfixed)
    ll(notfixed) = l;
end
l = ll;
uu = oldub;
if ~isempty(notfixed)
    uu(notfixed) = u;
end
u = uu;
infeasible = 0;

% Now find bounds from equalities
candidates = find(sum(Aeq | Aeq,2)==1);
if length(candidates)>0
    for i = candidates(:)'
        j = find(Aeq(i,:));
        new_bound = beq(i)/Aeq(i,j);
        if any(u(j) < new_bound) | any(l(j)> new_bound)
            infeasible = 1;
            return
        else
            u(j)=new_bound;
            l(j)=new_bound;
        end
    end
end

trivially_unbounded_below = find(trivially_unbounded_below);
trivially_unbounded_above = find(trivially_unbounded_above);
lb(trivially_unbounded_below(reset_to_minus_inf)) = -inf;
ub(trivially_unbounded_above(reset_to_inf)) = inf;
Matrices.lb = lb;
Matrices.ub = ub;
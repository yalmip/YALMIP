function p = remove_redundant(p);

if isempty(p.F_struc)
    return
end

b = p.F_struc(1+p.K.f:p.K.l+p.K.f,1);
A = -p.F_struc(1+p.K.f:p.K.l+p.K.f,2:end);

redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <-1e-2));

if length(redundant)>1
    p.InequalityConstraintState(redundant) = inf;
end

if p.options.bmibnb.lpreduce
    b = p.lpcuts(:,1);
    A = -p.lpcuts(:,2:end);
    redundant = find(((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <-1e-2));
    if length(redundant)>1
        p.lpcuts(redundant,:) = [];
        p.cutState(redundant) = [];
    end
end

if any(p.K.f)
    b = p.F_struc(1:p.K.f,1);
    A = -p.F_struc(1:p.K.f,2:end);
    s1 = ((A>0).*A*(p.ub-p.lb) - (b-A*p.lb) <1e-6);
    s2 = ((-A>0).*(-A)*(p.ub-p.lb) - ((-b)-(-A)*p.lb) <1e-6);
    redundant = find(s1 & s2);
    if length(redundant)>1
        p.EqualityConstraintState(redundant) = inf;
    end
end
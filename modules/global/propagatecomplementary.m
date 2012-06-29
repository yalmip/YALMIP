function p = propagatecomplementary(p)
LU = [p.lb p.ub];
complementary = find(p.lb==0 & p.ub==0 & p.variabletype'==1);
if ~isempty(complementary)
    for i = 1:length(complementary)
        index = find(p.bilinears(:,1) == complementary(i));
        x = p.bilinears(index,2);
        y = p.bilinears(index,3);
        if p.lb(x)>1e-4 | p.ub(x) < -1e-4 & p.lb(y)+p.ub(y)~=0
            if p.lb(y)>1e-4 | p.ub(y)<-1e-4
                p.feasible = 0;
            end
            p.ub(y)=0;
            p.lb(y)=0;
        elseif (p.lb(y)>1e-4 | p.ub(y) < -1e-4) & p.lb(x)+p.ub(x)~=0
            if p.lb(x)>1e-4 | p.ub(x)<-1e-4
                p.feasible = 0;
            end
            p.ub(x)=0;
            p.lb(x)=0;
        end
    end
end
complementary = find(p.lb==p.ub & p.variabletype'==1);
if ~isempty(complementary)
    for i = 1:length(complementary)
        index = find(p.bilinears(:,1) == complementary(i));
        x = p.bilinears(index,2);
        y = p.bilinears(index,3);
        k = p.lb(complementary(i));
        % x*y == k, case x and y negative, k positive
        if p.ub(x) < 0 & k >= 0
            p.lb(y) = max([p.lb(y) k/p.ub(x)]);
            p.ub(y) = min([p.ub(y) 0]);
        end
        if p.lb(y) > p.ub(y);p.feasible = 0;return;end
        if p.ub(y) < 0 & k >= 0
            p.lb(x) = max([p.lb(x) k/p.ub(y)]);
            p.ub(x) = min([p.ub(x) 0]);
        end
        if p.lb(x) > p.ub(x);p.feasible = 0;return;end
    end
end
if ~isequal(LU,[p.lb p.ub])
    p.changedbounds = 1;
end
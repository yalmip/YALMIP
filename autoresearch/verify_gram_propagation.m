function verify_gram_propagation
% Exhaustive feasible Gram assignments with randomized partial fixings.
rng(42,'twister');
n = 5;
m = n*(n-1)/2;
pair = zeros(n);
v = 0;
for j = 2:n
    for i = 1:j-1
        v = v+1;
        pair(i,j) = v;
        pair(j,i) = v;
    end
end
checked = 0;
for equalValue = 0:1
    for oppositeValue = [-1,-0.5]
        gram = {struct('pair',pair,'equalValue',equalValue, ...
            'saturated',true)};
        for pattern = 0:2^m-1
            z = (dec2bin(pattern,m)-'0')';
            X = eye(n);
            for j = 2:n
                for i = 1:j-1
                    if z(pair(i,j)) == equalValue
                        X(i,j) = 1;
                    else
                        X(i,j) = oppositeValue;
                    end
                    X(j,i) = X(i,j);
                end
            end
            if min(eig(X)) < -1e-10
                continue
            end
            for sample = 1:64
                fixed = rand(m,1) < 0.5;
                p = struct('lb',zeros(m,1),'ub',ones(m,1),'feasible',1);
                p.lb(fixed) = z(fixed);
                p.ub(fixed) = z(fixed);
                q = propagateGramEqualities(p,gram);
                assert(q.feasible && all(q.lb <= z) && all(q.ub >= z));
                checked = checked+1;
            end
        end
    end
end

gram = {struct('pair',pair,'equalValue',1,'saturated',true)};
p = struct('lb',zeros(m,1),'ub',ones(m,1),'feasible',1);
p.lb(pair(1,2)) = 1; p.ub(pair(1,2)) = 1;
p.lb(pair(1,3)) = 0; p.ub(pair(1,3)) = 0;
p.lb(pair(2,3)) = 1; p.ub(pair(2,3)) = 1;
q = propagateGramEqualities(p,gram);
assert(~q.feasible);

gram{1}.saturated = false;
p = struct('lb',zeros(m,1),'ub',ones(m,1),'feasible',1);
p.lb(1) = 1; p.ub(1) = 1;
[q,changed] = propagateGramEqualities(p,gram);
assert(~changed && isequal(q.lb,p.lb) && isequal(q.ub,p.ub));
fprintf('Verified %d feasible partial assignments.\n',checked);
end

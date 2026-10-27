function verify_gram_triangle_cuts
% Check structural detection and exhaustive validity of discrete Gram cuts.
n = 25;
nvars = n*(n-1)/2;
for equalValue = 0:1
    p.K = struct('f',0,'l',0,'q',[],'e',0,'p',[],'s',n);
    p.isbinary = true(1,nvars);
    p.F_struc = sparse(n*n,nvars+1);
    p.F_struc(1:n+1:n*n,1) = 1;
    pair = zeros(n);
    variable = 0;
    for j = 2:n
        for i = 1:j-1
            variable = variable+1;
            pair(i,j) = variable;
            pair(j,i) = variable;
            rows = [i+(j-1)*n,j+(i-1)*n];
            if equalValue == 1
                p.F_struc(rows,1) = -0.5;
                p.F_struc(rows,1+variable) = 1.5;
            else
                p.F_struc(rows,1) = 1;
                p.F_struc(rows,1+variable) = -1.5;
            end
        end
    end
    gram = initGramTriangleSeparator(p);
    assert(numel(gram) == 1 && gram{1}.equalValue == equalValue);

    % Enumerate every point on a five-vertex subgraph (ten binaries).
    small = pair(1:5,1:5);
    used = unique(nonzeros(small));
    [~,location] = ismember(small(small>0),used);
    small(small>0) = location;
    candidate = 0.5*ones(10,1);
    if equalValue == 1
        candidate([small(1,2),small(2,3),small(1,3)]) = [0.9,0.9,0.3];
    else
        candidate([small(1,2),small(2,3),small(1,3)]) = [0.1,0.1,0.7];
    end
    cuts = separateGramTriangles({struct('pair',small,'equalValue',equalValue)},candidate,10);
    assert(~isempty(cuts) && all(cuts*[1;candidate] < -1e-6));
    feasibleCount = 0;
    for pattern = 0:2^10-1
        z = (dec2bin(pattern,10)-'0')';
        X = eye(5);
        for j = 2:5
            for i = 1:j-1
                if z(small(i,j)) == equalValue
                    value = 1;
                else
                    value = -0.5;
                end
                X(i,j) = value;
                X(j,i) = value;
            end
        end
        if min(eig(X)) >= -1e-10
            assert(all(cuts*[1;z] >= -1e-10));
            feasibleCount = feasibleCount+1;
        end
    end
    fprintf('Equal value %d: %d cuts valid over %d feasible patterns.\n', ...
        equalValue,size(cuts,1),feasibleCount);
end
end

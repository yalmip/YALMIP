function verify_gram_turan
rng(42,'twister');
n = 40;
K = 3;
d = 1;
b = -1/(K-1);
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
g = {struct('pair',pair,'equalValue',1,'saturated',true, ...
    'diagonal',d,'otherEntry',b,'exactPattern',true)};
% A regular 40-simplex is PSD, and each pair takes a fractional Gram value.
offDiagonal = -d/(n-1);
X = eye(n)*d + (ones(n)-eye(n))*offDiagonal;
assert(min(eig(X)) > -1e-9);
x = ones(m,1)*(offDiagonal-b)/(d-b);
[cliqueCuts,cliqueViolation] = separateGramCliques(g,x,m,ones(m,1));
assert(isempty(cliqueCuts) && cliqueViolation < 1e-9);
[cuts,violation] = separateGramTuran(g,x,m,ones(m,1));
assert(~isempty(cuts) && violation > 0.3);
assert(min(cuts*[1;x]) < -0.3);
for sample = 1:500
    labels = randi(K,n,1);
    z = zeros(m,1);
    for j = 2:n
        for i = 1:j-1
            z(pair(i,j)) = labels(i) == labels(j);
        end
    end
    assert(min(cuts*[1;z]) >= -1e-9);
end
g{1}.equalValue = 0;
[oppositeCuts,oppositeViolation] = separateGramTuran(g,1-x,m,ones(m,1));
assert(~isempty(oppositeCuts) && oppositeViolation > 0.3);
for sample = 1:500
    labels = randi(K,n,1);
    z = ones(m,1);
    for j = 2:n
        for i = 1:j-1
            z(pair(i,j)) = labels(i) ~= labels(j);
        end
    end
    assert(min(oppositeCuts*[1;z]) >= -1e-9);
end
% A single variable may be reused on many edges. The coefficient must then
% count every edge, or the r-vertex bound could exclude the all-equal PSD
% assignment.
reused = g;
reused{1}.equalValue = 1;
reused{1}.pair = double(pair > 0);
[reusedCuts,~] = separateGramTuran(reused,0,1,1);
assert(~isempty(reusedCuts));
assert(all(reusedCuts*[1;1] >= -1e-9));
% SDPA decimals may put the intended equal entry just below the diagonal
% and the class limit just below an integer. The 3x3 witness and the
% all-other eigenvalue still certify the partition inequalities.
rounded = struct;
rounded.K = struct('f',0,'l',0,'q',[],'e',0,'p',[],'s',n);
rounded.isbinary = true(1,m);
rounded.F_struc = sparse(n*n,m+1);
rounded.F_struc(1:n+1:n*n,1) = d;
roundedOther = -0.166666667;
roundedEqual = d-1e-9;
for j = 2:n
    for i = 1:j-1
        rows = [i+(j-1)*n,j+(i-1)*n];
        rounded.F_struc(rows,1) = roundedOther;
        rounded.F_struc(rows,1+pair(i,j)) = roundedEqual-roundedOther;
    end
end
detected = initGramTriangleSeparator(rounded);
assert(numel(detected) == 1 && ~detected{1}.saturated && ...
    detected{1}.exactPattern);
[roundedClique,~] = separateGramCliques(detected,zeros(m,1),m,ones(m,1));
[roundedCount,~] = separateGramTuran(detected,zeros(m,1),m,ones(m,1));
assert(~isempty(roundedClique) && ~isempty(roundedCount));
roundedCuts = [roundedClique;roundedCount];
for sample = 1:500
    labels = randi(7,n,1);
    z = zeros(m,1);
    for j = 2:n
        for i = 1:j-1
            z(pair(i,j)) = labels(i) == labels(j);
        end
    end
    assert(min(roundedCuts*[1;z]) >= -1e-9);
end
fprintf('Verified PSD witness beyond clique cuts and 1,000 K-class assignments.\n');
fprintf('Verified rounded near-saturation cuts on 500 seven-class assignments.\n');
end

function test_robust_8

n = 3;
r = 4;
m = 4;

x = sdpvar(n,1);
w = sdpvar(m,1);
sdpvar t

A = randn(r,n);
b = randn(r,1);
C = randn(r,m);

A1 = randn(2*r,m/2);
A2 = randn(2*r,m/2);
b1 = rand(2*r,1)*5;
b2 = rand(2*r,1)*5;

[xc,R] = chebyball(A1*w(1:m/2) <= b1)
while isinf(R) | R==0
    A1 = randn(2*r,m/2);
    b1 = rand(2*r,1)*5;
    [xc,R] = chebyball(A1*w(1:m/2) <= b1)
end
[xc,R] = chebyball(A2*w(m/2+1:m) <= b2)
while isinf(R) | R==0
    A2 = randn(2*r,m/2);
    b2 = rand(2*r,1)*5;
    [xc,R] = chebyball(A2*w(m/2+1:end) <= b2)
end

% Exact automatic
objective = norm(A*x + b + C*w,1);
model = [uncertain(w), A1*w(1:m/2)<=b1, A2*w(m/2+1:m)<=b2];
[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','none'));
solvesdp(F1,h)
o1 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','affine'));
solvesdp(F1,h)
o2 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','projection','robust.lplp','enumeration'));
solvesdp(F1,h)
o3 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','projection','robust.lplp','duality'));
solvesdp(F1,h)
o4 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','enumeration'));
solvesdp(F1,h)
o5 = double(t);


mbg_asserttolequal(o4-o5,0, 1e-4);
mbg_asserttolequal(o3-o5,0, 1e-4);
mbg_asserttolequal(double(o2<o1),1, 1e-4);
mbg_asserttolequal(double(o3<o1),1, 1e-4);


% Exact automatic
model = [uncertain(w), norm(w,1) <= 1];
[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','none'));
solvesdp(F1,h)
o1 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','affine'));
solvesdp(F1,h)
o2 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','projection','robust.lplp','enumeration'));
solvesdp(F1,h)
o3 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','projection','robust.lplp','duality'));
solvesdp(F1,h)
o4 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','enumeration'));
solvesdp(F1,h)
o5 = double(t);


mbg_asserttolequal(o4-o5,0, 1e-5);
mbg_asserttolequal(o3-o5,0, 1e-5);
mbg_asserttolequal(double(o2<o1),1, 1e-5);
mbg_asserttolequal(double(o3<o1),1, 1e-5);


% Exact automatic
model = [uncertain(w), norm(w,1)+norm(w,1) <= 1];
[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','none'));
solvesdp(F1,h)
o1 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','affine'));
solvesdp(F1,h)
o2 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','projection','robust.lplp','enumeration'));
solvesdp(F1,h)
o3 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','projection','robust.lplp','duality'));
solvesdp(F1,h)
o4 = double(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','enumeration'));
solvesdp(F1,h)
o5 = double(t);

mbg_asserttolequal(o4-o5,0, 1e-5);
mbg_asserttolequal(o3-o5,0, 1e-5);
mbg_asserttolequal(double(o2<o1),1, 1e-5);
mbg_asserttolequal(double(o3<o1),1, 1e-5);

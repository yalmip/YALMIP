function test_robust_6

% Exact automatic
yalmip('clear')
n = 1;
r = 1;
m = 1;

x = sdpvar(n,1);
w = sdpvar(m,1);
sdpvar t

model = [uncertain(w), norm(w,2) <= 1];
objective = norm(2*x + 3 + 4*w,1);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','none'));
solvesdp(F1,h)
o1 = double(h);
[F1,h] = robustify([model,objective <= t],t,sdpsettings('robust.auxreduce','affine'));
solvesdp(F1,h)
o2 = double(h);
mbg_asserttolequal(o1-o2,0, 1e-5);
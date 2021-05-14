function tests = test_robust_8
tests = functiontests(localfunctions);

function test1(testCase)

ops = sdpsettings('verbose',0);

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

[xc,R] = chebyball(A1*w(1:m/2) <= b1);
while isinf(R) | R==0
    A1 = randn(2*r,m/2);
    b1 = rand(2*r,1)*5;
    [xc,R] = chebyball(A1*w(1:m/2) <= b1);
end
[xc,R] = chebyball(A2*w(m/2+1:m) <= b2);
while isinf(R) | R==0
    A2 = randn(2*r,m/2);
    b2 = rand(2*r,1)*5;
    [xc,R] = chebyball(A2*w(m/2+1:end) <= b2);
end

% Exact automatic
objective = norm(A*x + b + C*w,1);
model = [uncertain(w), A1*w(1:m/2)<=b1, A2*w(m/2+1:m)<=b2];
[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','none'));
optimize(F1,h,ops) ;
o1 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','affine'));
optimize(F1,h,ops) ;
o2 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','projection','robust.lplp','enumeration'));
optimize(F1,h,ops) ;
o3 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','projection','robust.lplp','duality'));
optimize(F1,h,ops) ;
o4 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','enumeration'));
optimize(F1,h,ops) ;
o5 = value(t);


testCase.assertTrue(abs(o4-o5) <= 1e-4);
testCase.assertTrue(abs(o3-o5) <= 1e-4);
testCase.assertTrue(o2<o1);
testCase.assertTrue(o3<o1);


% Exact automatic
model = [uncertain(w), norm(w,1) <= 1];
[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','none'));
optimize(F1,h,ops) ;
o1 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','affine'));
optimize(F1,h,ops) ;
o2 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','projection','robust.lplp','enumeration'));
optimize(F1,h,ops) ;
o3 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','projection','robust.lplp','duality'));
optimize(F1,h,ops) ;
o4 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','enumeration'));
optimize(F1,h,ops) ;
o5 = value(t);


testCase.assertTrue(abs(o4-o5) <= 1e-5);
testCase.assertTrue(abs(o3-o5) <= 1e-5);
testCase.assertTrue(o2<o1);
testCase.assertTrue(o3<o1);


% Exact automatic
model = [uncertain(w), norm(w,1)+norm(w,1) <= 1];
[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','none'));
optimize(F1,h,ops) ;
o1 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','affine'));
optimize(F1,h,ops) ;
o2 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','projection','robust.lplp','enumeration'));
optimize(F1,h,ops) ;
o3 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','projection','robust.lplp','duality'));
optimize(F1,h,ops) ;
o4 = value(t);

[F1,h] = robustify([model,objective <= t],t,sdpsettings('verbose',0,'robust.auxreduce','enumeration'));
optimize(F1,h,ops) ;
o5 = value(t);

testCase.assertTrue(abs(o4-o5) <= 1e-5);
testCase.assertTrue(abs(o3-o5) <= 1e-5);
testCase.assertTrue(o2<o1);
testCase.assertTrue(o3<o1);

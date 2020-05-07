function tests = test_robust_7
tests = functiontests(localfunctions);

function test1(dummy)
ops = sdpsettings;
yalmip('clear')
A = [2.938 -0.7345 0.25;4 0 0;0 1 0];
B = [0.25;0;0];
C = [-0.2072 0.04141 0.07256];
E = [0.0625;0;0];

N = 10;
U = sdpvar(N,1);
W = sdpvar(N,1);
x = sdpvar(3,1);

Y = [];
xk = x;
for k = 1:N
 xk = A*xk + B*U(k)+E*W(k);
 Y = [Y;C*xk];
end

F = [Y <= 1, -1 <= U <= 1];
objective = norm(Y-1,1) + norm(U,1)*0.01;

G = [-1 <= W <= 1]

% future controls are functions (policies) of future states, i.e. they
% depend on past disturbances w. 
% We can hack YALMIP to eliminate/adjust 
for i = 2:N
    yalmip('setdependenceUser',U(i),W(1:i-1))
end

% since we are comparing with manually constructed affine feedback, we must
% force YALMIP to use this. Note that it will be used on the real auxiliary
% variables too, i.e. sub-optimal
[Frobust,h] = robustify(F + G,objective,sdpsettings('robust.auxreduce','affine'),W);

Ns = 3;
xk1 = [0;0;0];
ww = (-1+2*rand(Ns,1))

for i = 1:Ns
    optimize([Frobust, x == xk1(:,end)],h,ops);    
    xk1 = [xk1 A*xk1(:,end) + B*value(U(1)) + E*ww(i)];
end
xk1

V = sdpvar(N,1);
L = sdpvar(N,N,'full').*(tril(ones(N))-eye(N));

U = L*W + V;

Y = [];
xk = x;
for k = 1:N
 xk = A*xk + B*U(k)+E*W(k);
 Y = [Y;C*xk];
end

F = [Y <= 1, -1 <= U <= 1];
objective = norm(Y-1,1) + norm(U,1)*0.01;

[Frobust,h] = robustify([F, G],objective,sdpsettings('robust.auxreduce','affine'),W);

xk2 = [0;0;0];
ops = sdpsettings;
for i = 1:Ns
    optimize([Frobust, x == xk2(:,end)],h,ops);
    xk2 = [xk2 A*xk2(:,end) + B*value(U(1)) + E*ww(i)];
end
norm(xk1-xk2)
assert(norm(xk1-xk2) <= 1e-3)




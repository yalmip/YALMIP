n = 10;
m = 10;

Q = randn(n,n);Q = Q*Q';
c = randn(n,1);
d = randn(m,1);
A = randn(15,n);
b = rand(15,1)*20*n;
E = randn(15,m);
S1 = sprand(n+m,1,0.5);
S2 = sprand(n+m,1,0.5);

H = randn(m,m);H = H*H';
e = randn(m,1);
f = randn(n,1);
F = randn(10,m);
h = rand(10,1)*2*m;
G = randn(10,n);


x = sdpvar(n,1);
z = sdpvar(m,1);
lambda = sdpvar(length(h),1);
nu = sdpvar(1)
slack = h + G*x - F*z;
KKT = [H*z + e + F'*lambda+S1(n+1:end)*nu == 0,
    F*z <=  h + G*x,
    S1'*[x;z] == 0,
    lambda >= 0];

for i = 1:length(h)
    KKT = [KKT, ((lambda(i)==0) | (slack(i) == 0))];
end
KKT = [KKT, lambda <= 1000, -1000 <= [x;z] <= 1000];

sol = solvesdp([KKT, A*x <= b + E*z,S2'*[x;z] == 0], 0.5*x'*Q*x + c'*x + d'*z);
%t = sdpvar(1)
%sol = solvesdp([KKT, A*x <= b + E*z,S2'*[x;z] == 0, 0.5*x'*Q*x + c'*x + d'*z<t],t);
if sol.problem == 0
    x1 = double(x);
    z1 = double(z);

    CO = [A*x <= b+E*z,S2'*[x;z] == 0];%,0.5*x'*Q*x + c'*x + d'*z<t];
    CI = [F*z <= h + G*x,S1'*[x;z] == 0];
    OO = 0.5*x'*Q*x + c'*x + d'*z;
    OI = 0.5*z'*H*z + e'*z + f'*x;
    [sol,xsol,zsol,info] = solvebilevel2(CO,OO,CI,OI,z);
    xsol-x1
    zsol-z1
end


% Lp from Bard, with integer variables
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2-4*y3;
CO = [[x1 x2]>0, integer([x1 x2])];
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] > 0,-y1+y2+y3 < 10, 2*x1-y1+2*y2-0.5*y3 < 10, 2*x2+2*y1-y2-0.5*y3 < 9.7]
sol = solvebilevel(CO,OO,CI,OI,[y1 y2 y3]);
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3]);
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.rootcuts',20));



% Lp from bard, no y3 in outer problem
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2;
CO = [x1 x2]>0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] > 0,-y1+y2+y3 < 1, 2*x1-y1+2*y2-0.5*y3 < 1, 2*x2+2*y1-y2-0.5*y3 < 1]
solvebilevel(CO,OO,CI,OI,[y1 y2 y3]);
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.rootcuts',20));


% Lp from bard, no y3 in outer problem
sdpvar x1 x2 y1 y2 y3
OO = -8*x1-4*x2+4*y1-40*y2+x1^2+x2^2;
CO = [x1 x2]>0;
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] > 0,-y1+y2+y3 < 1, 2*x1-y1+2*y2-0.5*y3 < 1, 2*x2+2*y1-y2-0.5*y3 < 1]
[sol,xsol,zsol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3]);
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',0));
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.rootcuts',20));


% Lp from bard, with added quadratic term handled by SDP
sdpvar x1 x2 y1 y2 y3 t1 t2
OO = -8*x1-4*x2+4*y1-40*y2+t1+t2;
CO = [[x1 x2]>0,[t1 x1;x1 1]>0,[t2 x2;x2 1]>0];
OI = x1+2*x2+y1+y2+2*y3;
CI = [[y1 y2 y3] > 0,-y1+y2+y3 < 1, 2*x1-y1+2*y2-0.5*y3 < 1, 2*x2+2*y1-y2-0.5*y3 < 1]
[sol,xsol,zsol,info] = solvebilevel(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',1));
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.solvefrp',0));
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.rootcuts',20));

n = 10;
m = 10;
r = 40;
x = sdpvar(n,1);
y = sdpvar(m,1);
b1 = rand(r,1)*r*3+1;
b2 = rand(r,1)*r*3+1;
c1 = randn(n,1);
c2 = randn(m,1);
d1 = randn(n,1);
d2 = randn(m,1);
A1 = randn(r,n);
B1 = randn(r,m);
A2 = randn(r,n);
B2 = randn(r,m);
Q2 = randn(n+m);Q2 = Q2*Q2';

OO = c1'*x+d1'*y;
CO = [A1 B1]*[x;y] < b1;
OI = c2'*x+d2'*y+[x;y]'*Q2*[x;y];
CI = [A2 B2]*[x;y] < b2;
solvebilevel(CO,OO,CI,OI,y);
solvebilevel(CO,OO,CI,OI,y,sdpsettings('bilevel.solvefrp',0));
lambda = sdpvar(length(b2),1);
s = sdpvar(length(b2),1);
solvesdp([CO,CI,200>lambda>0,s == b2-[A2 B2]*[x;y],lambda.*s==0,d2+B2'*lambda==0],OO,sdpsettings('solver','bmibnb','bmibnb.upper','snopt','debug',1,'bmibnb.maxiter',200))
solvesdp([CO,CI,200>lambda>0,s == b2-[A2 B2]*[x;y],d2-B2'*lambda==0],OO,sdpsettings('solver','bmibnb','bmibnb.upper','snopt','debug',1,'bmibnb.maxiter',200))



lambda = sdpvar(length(b2),1);
s = sdpvar(length(b2),1);
solvesdp([CO,CI,lambda>0,s == b2-[A2 B2]*[x;y], lambda.*s==0,d2-B2'*lambda==0],OO,sdpsettings('solver','bmibnb'))



% Bard falk
sdpvar x1 x2 y1 y2 y3
OO = -8*x1 - 4*x2 + 4*y1 - 40*y2 + 4*y3;
OI = x1 + 2*x2 + y1 + y2 + 2*y3;
CI = [-y1 + y2 + y3 < 1,2*x1 - y1 + 2*y2 - 0.5*y3 < 1,2*x2 + 2*y1 - y2 - 0.5*y3 < 1,[y1 y2 y3]>0,[x1 x2]>0];
CO = [[y1 y2 y3]>0,[x1 x2]>0];
solvebilevel(CO,OO,CI,OI,[y1 y2 y3]);
% This one fails with CLP!
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y1 y2 y3],sdpsettings('bilevel.rootcuts',20,'bilevel.solvefrp',0));



% Audet Hafddad Savard
sdpvar x y 
OO = x+2*y
CO = x>0
OI = -4*x-y
CI = [x+y>8,3*x-2*y>-6, -3*x-4*y>-48,y>0]
solvebilevel(CO,OO,CI,OI,[y]);
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y],sdpsettings('bilevel.rootcuts',inf));


sdpvar x y
OO = -x-10*y;
CO = [x>0,y>0];
OI = y;
CI = [-25*x+20*y<30,x+2*y<10,2*x-y < 15,2*x+10*y>15,x>0,y>0];
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y]);
[sol,xsol,zsol,info] = solvebileveldisjunction(CO,OO,CI,OI,[y],sdpsettings('bilevel.rootcuts',inf));




A = randn(25,3);
x = [1;2;3];

e = sprand(25,1,0.2);
y = A*x+e


theta = sdpvar(3,1);

solvesdp(theta>0,(y-A*theta)'*(y-A*theta))

delta = sdpvar(25,1);
CI = theta>0
CO = [-1 < delta<1,theta > 0]
OI = (y+delta-A*theta)'*(y+delta-A*theta);
OO = 0.1*norm(delta,1)+(y-A*theta)'*(y-A*theta);
[sol,xsol,zsol,info] = solvebilevel2(CO,OO,CI,OI,theta);

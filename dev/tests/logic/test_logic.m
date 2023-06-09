function tests = test_logic_constraint_eq_1
tests = functiontests(localfunctions);

function test1(testCase)
binvar a b
intvar x
ops=sdpsettings('verbose',0);
F = [ a == (b & (x<=4.5)), -15 <= x <= 15];
sol = optimize(F,a^2,ops);
testCase.assertTrue(sol.problem == 0);
% Since a is 0, b has to be false or constraint violated
testCase.assertTrue(value(a)==0 && (value(b)==0 || value(x)>=5));

F = true(or(0,a));
sol = optimize(F,a,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(value(a)==1);

F = true(and(0,a));
sol = optimize(F,a,ops);
testCase.assertTrue(sol.problem == 1);

F = true(xor(0,a));
sol = optimize(F,a,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(value(a)==1);

F = true(xor(1,a));
sol = optimize(F,-a,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(value(a)==0);

F = [ a == (b & (x<=4.5)), -15 <= x <= 15];
sol = optimize(F,x^2 + 1-b,ops);
testCase.assertTrue(value(a)==1 && (value(b)==1 || value(x)==0));

dX = binvar(1);
dY = binvar(1);
sys = [iff(dX,a), iff(dY,[b & (x<=4.5)]), dX == dY];
F = [sys, -15 <= x <= 15];
sol = optimize(F,x^2 + 1-b,ops);
testCase.assertTrue(value(a)==1 && (value(b)==1 || value(x)==0));

a = binvar(3,2);
b = binvar(3,2);
optimize([b(1)==0, true(or(a,b))],sum(sum(a-b)),ops);
testCase.assertTrue(value(sum(sum(a-b)))==-4);

optimize([b(1)==1, true(and(a,b))],sum(sum(a-b)),ops);
testCase.assertTrue(value(sum(sum(a-b)))==0);

optimize([b(1)==1, true(xor(a,b))],sum(sum(-a-b)),ops);
testCase.assertTrue(value(sum(sum(-a-b)))==-6);

optimize([true(a & (~b))],sum(sum(a-b)),ops);
testCase.assertTrue(value(sum(sum(a-b)))==6);

function test2(testCase)
yalmip('clear')
ops=sdpsettings('verbose',0);

binvar x y
optimize(iff(x,y <= 0.5),2*x+y,ops);
testCase.assertTrue(abs(value(x)) <= 1e-4);
testCase.assertTrue(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y <= 0.5),x+2*y,ops);
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(abs(value(y)) <= 1e-4);

optimize(iff(x,y == 1),x+y,ops);
testCase.assertTrue(abs(value(x))<=1e-4)
testCase.assertTrue(abs(value(y))<=1e-4)

optimize(iff(x,y == 1),-2*x+y,ops);
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y == 1),x-2*y,ops);
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(abs(value(y)-1)<= 1e-4);

optimize(iff(x,y == 0),x,ops);
testCase.assertTrue(abs(value(x))<=1e-4)
testCase.assertTrue(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y == 0),y,ops);
testCase.assertTrue(abs(value(x)-1) <= 1e-4);
testCase.assertTrue(abs(value(y))<=1e-4)

optimize(iff(x,y == 0),2*x+y,ops);
testCase.assertTrue(abs(value(x))<=1e-4)
testCase.assertTrue(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y == 0),x+2*y,ops);
testCase.assertTrue(abs(value(x)-1)<=1e-4)
testCase.assertTrue(abs(value(y))<=1e-4)

% Polytopic i.f.f binary equality
sdpvar x
binvar y
optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],y,ops);
testCase.assertTrue(abs(value(x))>.0999);
testCase.assertTrue(abs(value(y)) <= 1e-3);
optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x^2,ops);
testCase.assertTrue(abs(value(x))<=1e-3)
testCase.assertTrue(abs(value(y)-1)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],200*x^2+y,ops);
testCase.assertTrue(abs(value(x))<=1e-3)
testCase.assertTrue(abs(value(y)-1)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],-x,ops);
testCase.assertTrue(abs(value(x)-3)<=1e-3)
testCase.assertTrue(abs(value(y))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x,ops);
testCase.assertTrue(abs(value(x)--3)<=1e-3)
testCase.assertTrue(abs(value(y))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x-y,ops);
testCase.assertTrue(abs(value(x)--3)<=1e-3)
testCase.assertTrue(abs(value(y))<=1e-3)


% Polytopic x equality iff binary
sdpvar x
binvar y z
optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3<=x<=3,-10<=z<=10],y,ops);
testCase.assertTrue(abs(value(y))<=1e-3)
testCase.assertTrue(abs(value(z))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-y,ops);
testCase.assertTrue(abs(value(y)-1)<=1e-3)
testCase.assertTrue(abs(value(z)-1)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x,ops);
testCase.assertTrue(abs(value(y))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x-y,ops);
testCase.assertTrue(abs(value(y))<=1e-3)
testCase.assertTrue(abs(value(x)--3)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x-5*y,ops);
testCase.assertTrue(abs(value(y)-1)<=1e-3)
testCase.assertTrue(abs(value(z)-1)<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-y,ops);
testCase.assertTrue(abs(value(y)-1)<=1e-3)
testCase.assertTrue(abs(value(z)-1)<=1e-3)
testCase.assertTrue(abs(value(x)-2)<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-z+x,ops);
testCase.assertTrue(abs(value(y))<=1e-3)
testCase.assertTrue(abs(value(z)-1)<=1e-3)
testCase.assertTrue(abs(value(x)--3)<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],z,ops);
testCase.assertTrue(abs(value(y))<=1e-3)
testCase.assertTrue(abs(value(z))<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],z+x,ops);
testCase.assertTrue(abs(value(y))<=1e-3)
testCase.assertTrue(abs(value(z))<=1e-3)
testCase.assertTrue(abs(value(x)--3)<=1e-3)

sdpvar z
optimize([iff([x <= 2, z >= 1],y == 1),-3 <= x <= 3,-10<=z<=10],z,ops);
testCase.assertTrue(abs(value(y))<=1e-3)
testCase.assertTrue(abs(value(z)--10)<=1e-3)

optimize([iff([x <= 2, z >= 1],y == 1),-3 <= x <= 3,-10<=z<=10],-z+x,ops);
testCase.assertTrue(abs(value(y)-1)<=1e-3)
testCase.assertTrue(abs(value(z)-10)<=1e-3)
testCase.assertTrue(abs(value(x)--3)<=1e-3)

function test3(testCase)
binvar d1 d2
ops=sdpsettings('verbose',0);

sol =  optimize(~(d1==d2),d2,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(d1 + d2)-1) <= 1e-3)

sol =  optimize((d1~=d2),d1+d2,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(d1 + d2)-1) <= 1e-3)

sol =  optimize(iff(d1,d2),d1+d2,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(d1 + d2)-0) <= 1e-3)

sol =  optimize(iff(d1,d2),d1,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(abs(value(d2)-0) <= 1e-3)

function test4(testCase)
yalmip('clear')
ops=sdpsettings('verbose',0);

binvar x y
optimize([0 <= [x,y] <= 1, x+y~=1],(x+y-1)^2-10*x,ops);
testCase.assertTrue(abs(value(x)-1) <= 1e-4)
testCase.assertTrue(abs(value(y)-1) <= 1e-4)

binvar x y
optimize([0 <= [x,y] <= 1, x+y~=1],(x+y-1)^2,ops);
testCase.assertTrue(abs(value((x+y-1)^2)-1) <= 1e-4)

binvar x1 y1 x2 y2
optimize([0 <= [x1,y1,x2,y2] <= 1, [x1;x2]'+[y1;y2]'~=1],(x1+y1-1)^2-10*x1+(x2+y2-1)^2-10*x2,ops);
testCase.assertTrue(abs(value(x1)-1) <= 1e-4)
testCase.assertTrue(abs(value(y1)-1) <= 1e-4)

function test5(testCase)
ops=sdpsettings('verbose',0);

x = sdpvar(4,1);
sol = optimize((-10 <= x <= 10) + (nnz(x)==3),(x-2)'*(x-2),ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(norm(sort(value(x))-sort([2;2;2;0])) < 1e-3)

sol = optimize((-10 <= x <= 10) + (nnz(x)<=3),(x-2)'*(x-2),ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(norm(sort(value(x))-sort([2;2;2;0])) < 1e-3)

sdpvar x y z
sol = optimize((200>=[x y z] >= -20) + (nnz([x;y;z] == 3) >= 2),2*x+y+z,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(norm(value([x y z] - [-20 3 3])) < 1e-3)


function test7(testCase)
ops=sdpsettings('verbose',0);

binvar aa ba ca da ea fa % NOTE a b c ... generates error when run as function !
F = (true((aa & ba & ca) | (da & ea & fa)));
sol = optimize(F,[],ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue( all(value([aa ba ca])== [1 1 1])  | all(value([da ea fa])== [1 1 1]) )


function test8(testCase)
ops=sdpsettings('verbose',0);

randn('seed',12345);
rand('seed',12345);

A1 = randn(8,2);
b1 = rand(8,1)*2-A1*[3;3];
A2 = randn(8,2);
b2 = rand(8,1)*2-A2*[-3;3];
A3 = randn(8,2);
b3 = rand(8,1)*2-A3*[3;-3];
A4 = randn(8,2);
b4 = rand(8,1)*2-A4*[-3;-3];

binvar inp1 inp2 inp3 inp4
F = true(inp1 | inp2 | inp3 | inp4);
x = sdpvar(2,1);
F = F + (iff(inp1,A1*x <= b1));
F = F + (iff(inp2,A2*x <= b2));
F = F + (iff(inp3,A3*x <= b3));
F = F + (iff(inp4,A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = optimize(F,-x(2),ops);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(norm(value(x) - [1.89084694511033   4.32311245145436]') <= 1e-6)
testCase.assertTrue(norm(value([inp1 inp2 inp3 inp4]) - [0 0 0 1]) <= 1e-6)

F = true(inp1 | inp2 | inp3 | inp4);
F = F + (inp1 == (A1*x <= b1));
F = F + (inp2 == (A2*x <= b2));
F = F + (inp3 == (A3*x <= b3));
F = F + (inp4 == (A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = optimize(F,-x(2),ops);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(norm(value(x) - [1.89084694511033   4.32311245145436]') <= 1e-6)
testCase.assertTrue(norm(value([inp1 inp2 inp3 inp4]) - [0 0 0 1]) <= 1e-6)

F = true(inp1 | inp2 | inp3 | inp4);
F = F + (implies(inp1,A1*x <= b1));
F = F + (implies(inp2,A2*x <= b2));
F = F + (implies(inp3,A3*x <= b3));
F = F + (implies(inp4,A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = optimize(F,-x(2),ops);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(norm(value(x) - [1.89084694511033   4.32311245145436]') <= 1e-6)
testCase.assertTrue(norm(value([inp1 inp2 inp3 inp4]) - [0 0 0 1]) <= 1e-6)

F = ( (A1*x <= b1) | (A2*x <= b2) | (A3*x <= b3) | (A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = optimize(F,-x(2),ops);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(norm(value(x) - [1.89084694511033   4.32311245145436]') <= 1e-6)
testCase.assertTrue(norm(value([inp1 inp2 inp3 inp4]) - [0 0 0 1]) <= 1e-6)

x = sdpvar(2,1);
bounds(x,-100,100);
F = ( (A1*x <= b1) | (A2*x <= b2) | (A3*x <= b3) | (A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = optimize(F,-x(2),ops);
testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(norm(value(x) - [1.89084694511033   4.32311245145436]') <= 1e-6)
testCase.assertTrue(norm(value([inp1 inp2 inp3 inp4]) - [0 0 0 1]) <= 1e-6)

ii = sdpvar(1,1);
jj = sdpvar(1,1);
x = sdpvar(1,8);
p = [0 1 7 2 3 4 3 20];
optimize((-100 <= [x(:);ii;jj] <= 100) + (x == p)+(x([ii jj]) <= 3)+(ii~=jj),-ii-jj,ops);
testCase.assertTrue(sol.problem == 0);
testCase.assertTrue(fix(min(value([ii jj]))) == [5]);
testCase.assertTrue(fix(max(value([ii jj]))) == [7]);




function tests = test_logic_iff
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear')

binvar x y
optimize(iff(x,y <= 0.5),2*x+y)
assert(abs(value(x)) <= 1e-4);
assert(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y <= 0.5),x+2*y)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(y)) <= 1e-4);

optimize(iff(x,y == 1),x+y)
assert(abs(value(x))<=1e-4)
assert(abs(value(y))<=1e-4)

optimize(iff(x,y == 1),-2*x+y)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y == 1),x-2*y)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(y)-1)<= 1e-4);

optimize(iff(x,y == 0),x)
assert(abs(value(x))<=1e-4)
assert(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y == 0),y)
assert(abs(value(x)-1) <= 1e-4);
assert(abs(value(y))<=1e-4)

optimize(iff(x,y == 0),2*x+y)
assert(abs(value(x))<=1e-4)
assert(abs(value(y)-1)<= 1e-4);
optimize(iff(x,y == 0),x+2*y)
assert(abs(value(x)-1)<=1e-4)
assert(abs(value(y))<=1e-4)

% Polytopic i.f.f binary equality
sdpvar x
binvar y
optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],y,sdpsettings('debug',1))
assert(abs(value(x))>.0999);
assert(abs(value(y)) <= 1e-3);
optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x^2,sdpsettings('debug',1))
assert(abs(value(x))<=1e-3)
assert(abs(value(y)-1)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],200*x^2+y,sdpsettings('debug',1))
assert(abs(value(x))<=1e-3)
assert(abs(value(y)-1)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],-x,sdpsettings('debug',1))
assert(abs(value(x)-3)<=1e-3)
assert(abs(value(y))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x,sdpsettings('debug',1))
assert(abs(value(x)--3)<=1e-3)
assert(abs(value(y))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1],y == 1),-3<=x<=3],x-y,sdpsettings('debug',1))
assert(abs(value(x)--3)<=1e-3)
assert(abs(value(y))<=1e-3)


% Polytopic x equality iff binary
sdpvar x
binvar y z
optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3<=x<=3,-10<=z<=10],y,sdpsettings('debug',1))
assert(abs(value(y))<=1e-3)
assert(abs(value(z))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-y,sdpsettings('debug',1))
assert(abs(value(y)-1)<=1e-3)
assert(abs(value(z)-1)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x,sdpsettings('debug',1))
assert(abs(value(y))<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x-y,sdpsettings('debug',1))
assert(abs(value(y))<=1e-3)
assert(abs(value(x)--3)<=1e-3)

optimize([iff([-0.1 <= x <= 0.1, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],x-5*y,sdpsettings('debug',1))
assert(abs(value(y)-1)<=1e-3)
assert(abs(value(z)-1)<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-y,sdpsettings('debug',1))
assert(abs(value(y)-1)<=1e-3)
assert(abs(value(z)-1)<=1e-3)
assert(abs(value(x)-2)<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],-z+x,sdpsettings('debug',1))
assert(abs(value(y))<=1e-3)
assert(abs(value(z)-1)<=1e-3)
assert(abs(value(x)--3)<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],z,sdpsettings('debug',1))
assert(abs(value(y))<=1e-3)
assert(abs(value(z))<=1e-3)

optimize([iff([x == 2, z==1],y == 1),-3 <= x <= 3,-10<=z<=10],z+x,sdpsettings('debug',1))
assert(abs(value(y))<=1e-3)
assert(abs(value(z))<=1e-3)
assert(abs(value(x)--3)<=1e-3)

sdpvar z
optimize([iff([x <= 2, z >= 1],y == 1),-3 <= x <= 3,-10<=z<=10],z,sdpsettings('debug',1))
assert(abs(value(y))<=1e-3)
assert(abs(value(z)--10)<=1e-3)

optimize([iff([x <= 2, z >= 1],y == 1),-3 <= x <= 3,-10<=z<=10],-z+x,sdpsettings('debug',1))
assert(abs(value(y)-1)<=1e-3)
assert(abs(value(z)-10)<=1e-3)
assert(abs(value(x)--3)<=1e-3)


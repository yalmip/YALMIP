function tests = test_sos_book
tests = functiontests(localfunctions);

function test1(dummy)

k=2;
mu = sdpvar
sdpvar x y 
[h,c]=polynomial([x y],2*k-4);
g=y^2*(1-x^2)-(x^2+2*y-1)^2;
F=sos(y+mu-h*g);
sol = solvesos(F,mu,[],[mu;c]);

assert(sol.problem == 0);
assert(abs(value(mu) - 0.17759) <= 1e-3);

k=4;
sdpvar x y 
[h,c]=polynomial([x y],2*k-4);
g=y^2*(1-x^2)-(x^2+2*y-1)^2;
F=sos(y+mu-h*g);
solvesos(F,mu,[],[mu;c]);
assert(sol.problem == 0);
assert(abs(value(mu) - 0.0161) <= 1e-3)

sdpvar x y lambda
f=x^4+y^4+2*x^2*y^2-x^2+y^2;
l=1/sqrt(8)-y;
F=sos(l+lambda*f);
sol = solvesos(F,0,[],lambda);

assert(sol.problem==0)
assert(length(sosd(F))==6)

sdpvar y1 y2 y3
M=[1 y1 y2;
y1 y2 y3;
y2 y3 y3+y2-y1];
sol = optimize(M>=0,y1);
assert(sol.problem==0);
sol = optimize(M>=0,-y1);
assert(sol.problem==0);

function tests = test_sos_matrix_1
tests = functiontests(localfunctions);

function test1(dummy)

sdpvar x y
P = [1+x^2 -x+y+x^2;-x+y+x^2 2*x^2-2*x*y+y^2];
m = size(P,1);
v = monolist([x y],degree(P)/2);
Q = sdpvar(length(v)*m);
R = kron(eye(m),v)'*Q*kron(eye(m),v)-P;
s = coefficients(R(findelements(triu(R))),[x y]);
sol = optimize((Q >= 0) + (s==0));
diff = (clean(P - kron(eye(m),v)'*value(Q)*kron(eye(m),v),1e-6));

assert(sol.problem == 0)
assert(isequal(diff,[0 0;0 0]))

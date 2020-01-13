function tests = test_operator_optimizer13
tests = functiontests(localfunctions);

function test1(dummy)

% Testing nonlinearly parameterized SDP

A = sdpvar(3,3,'full');
P = sdpvar(3,3);

A0 = randn(3);A0 = -A0*A0';
C =  [A'*P+P*A <= -eye(3), P>=0];
M = optimizer(C,trace(P),sdpsettings('solver','+mosek'),A,P);
P0 = M{A0};
optimize([A0'*P+P*A0 <= -eye(3), P>=0],trace(P));
assert(abs(trace(P0)-trace(value(P))) <= 1e-3)

Q0 = randn(3);Q0 = Q0*Q0';
Q = sdpvar(3);
M = optimizer(C,trace(Q*P),sdpsettings('solver','+mosek'),{A,Q},P);
P0 = M{{A0,Q0}};
optimize([A0'*P+P*A0 <= -eye(3), P>=0],trace(Q0*P));
assert(abs(trace(Q0*P0)-trace(Q0*value(P))) <= 1e-3)

% Test a case where an SDP constraint boils down to a semidefinite constant
sdpvar x y
P = optimizer([[y*x y;y 1] >=0, [x 1;1 2]>=0],x,sdpsettings('solver','+mosek'),y,x);
assert(abs(P{0}-.5) <= 1e-4);
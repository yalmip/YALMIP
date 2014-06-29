function test_operator_optimizer13

% Testing nonlinearly parameterized SDP

A = sdpvar(3,3,'full');
P = sdpvar(3,3);

A0 = randn(3);A0 = -A0*A0';
C =  [A'*P+P*A <= -eye(3), P>=0];
M = optimizer(C,trace(P),sdpsettings('solver','+sdpt3'),A,P);
P0 = M{A0};
solvesdp([A0'*P+P*A0 <= -eye(3), P>=0],trace(P));
mbg_asserttolequal(trace(P0),trace(double(P)),1e-3)

Q0 = randn(3);Q0 = Q0*Q0';
Q = sdpvar(3);
M = optimizer(C,trace(Q*P),sdpsettings('solver','+sdpt3'),{A,Q},P);
P0 = M{{A0,Q0}};
solvesdp([A0'*P+P*A0 <= -eye(3), P>=0],trace(Q0*P));
mbg_asserttolequal(trace(Q0*P0),trace(Q0*double(P)),1e-3)

% Test a case where an SDP constraint boils down to a semidefinite constant
sdpvar x y
P = optimizer([[y*x y;y 1] >=0, [x 1;1 2]>=0],x,sdpsettings('solver','+sdpt3'),y,x);
mbg_asserttolequal(P{0},.5,1e-4);
function tests = test_mpt_oneshotmpqp
tests = functiontests(localfunctions);

function test1(dummy)

yalmip('clear')

% Model data
A = [2 -1;1 0];
B = [1;0];
C = [0.5 0.5];
nx = 2; % Number of states
nu = 1; % Number of inputs

% Prediction horizon
N = 5;
% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));
% Inputs u(k), ..., u(k+N) (last one not used) u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

F = ([]);
obj = 0;
for k = N-1:-1:1
    % Feasible region
    F = F + (-1 <= u{k}     <= 1);
    F = F + (-1 <= C*x{k}   <= 1);
    F = F + (-5 <= x{k}     <= 5);
    F = F + (-1 <= C*x{k+1} <= 1);
    F = F + (-5 <= x{k+1}   <= 5);
    % Dynamics
    F = F + (x{k+1} == A*x{k}+B*u{k});
    % Cost in value iteration
    obj = obj + x{k}'*x{k} + u{k}'*u{k};
end

[sol{k},diagnost{k},Uz{k},J{k},Optimizer{k}] = solvemp(F,obj,[],x{k},u{k});
assign(x{k},[1;0.5])

assert(diagnost{1}.problem == 0);
assert(abs(value(J{k})-3.82456140350877) <= 1e-5)
assert(abs(value(Optimizer{k})--0.99122807017544) <= 1e-5)
assign(x{k},[0;1.7])
assert(abs(value(Optimizer{k})-1) <= 1e-5)
assert(abs(value(J{k})-6.19300000000000) <= 1e-5)



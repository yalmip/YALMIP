function tests = test_operator_optimizer10
tests = functiontests(localfunctions);

function test1(dummy)
% Model data
A0 = [2 -1;1 0.207];
A = sdpvar(2,2,'full');
B = [1;.12];
B0 = B;
nx = 2; % Number of states
nu = 1; % Number of inputs

% MPC data
Q = eye(2);
R = 2;
N = 8;

u = sdpvar(repmat(nu,1,N),repmat(1,1,N));
x = sdpvar(repmat(nx,1,N+1),repmat(1,1,N+1));

constraints = [];
objective = 0;
for k = 1:N
    objective = objective + norm(Q*x{k},1) + norm(R*u{k},1);
    constraints = [constraints, x{k+1} == A*x{k} + B*u{k}];
    constraints = [constraints, -1 <= u{k}<= 1, -5<=x{k}<=5];
end

controller = optimizer(constraints, objective,sdpsettings('verbose',1,'solver','+bnb'),{A,x{1}},u{1});
x0 = [3;1]*.25;
m = 2;
all = {repmat(A0,1,m)+repmat(.5*[0 1;.3 1],1,m).*randn(2,2*m)*.02, repmat(x0,1,m)+randn(2,m)*0};
U = controller{all};
assert(isequal(size(U),[1 2]));

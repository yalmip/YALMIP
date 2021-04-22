function tests = test_mpt_lti_1
tests = functiontests(localfunctions);

function test1(dummy)

% Data
A = [2 -1;1 0];nx = 2;
B = [1;0];nu = 1;
C = [0.5 0.5];

% Prediction horizon
N = 3;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Value functions
J = cell(1,N);

% Initialize value function at stage N
J{N} = 0;
U{N} = 0;
k=N-1;

for k = N-1:-1:1
    
    % Feasible region       
    bounds(x{k},-5,5);
    bounds(u{k},-1,1);
    bounds(x{k+1},-5,5);
    
    F =     (-1 <= u{k}     <= 1);
    F = F + (-1 <= C*x{k}   <= 1);
    F = F + (-5 <= x{k}     <= 5);
    F = F + (-1 <= C*x{k+1} <= 1);
    F = F + (-5 <= x{k+1}   <= 5);  
  
    % LTI Dynamics
    F = F + (x{k+1} == A*x{k}+B*u{k});
    
    obj = norm([x{k};u{k}],1);
    
    % Compute value function for one step backwards
    [mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,obj + J{k+1},[],x{k},u{k});    
end

% MPT implementation
sysStruct.A= A;
sysStruct.B= B;
sysStruct.C= C;
sysStruct.D= [0];

%set constraints on output
sysStruct.ymin    =   -1;
sysStruct.ymax    =    1;

%set constraints on input
sysStruct.umin    =   -1;
sysStruct.umax    =   1;

sysStruct.xmin    =   [-5;-5];
sysStruct.xmax    =   [5;5];

probStruct.norm=1;
probStruct.Q=eye(2);
probStruct.R=1;
probStruct.N=N-1;
probStruct.P_N=zeros(2);
probStruct.subopt_lev=0;
probStruct.Tconstraint=0;

ctrl = mpt_control(sysStruct, probStruct, 'online');
Y = ctrl.toYALMIP();
Y.constraints = Y.constraints + [ -1 <= C*Y.variables.x(:, end) <= 1 ];
new = ctrl.fromYALMIP(Y).toExplicit();

fun1 = mpt_mpsol2pu(mpsol{1});
result = fun1.compare(new.optimizer, 'obj');
assert(result == 0);

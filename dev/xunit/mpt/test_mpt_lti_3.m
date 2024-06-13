function tests = test_mpt_lti_3
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
    
    obj =x{k}'*x{k} + u{k}'*u{k};
    
    % Compute value function for one step backwards
    [mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,obj + J{k+1},[],x{k},u{k});    
end

assert(sol{1}.problem == 0);

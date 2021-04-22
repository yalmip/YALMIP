function tests = test_mpt_pwa_1
tests = functiontests(localfunctions);

function test1(dummy)

% Data
A = [2 -1;1 0];nx = 2;
B = [1;0];nu = 1;
C = [0.5 0.5];

% Prediction horizon
N = 6;

% States x(k), ..., x(k+N)
x = sdpvar(repmat(nx,1,N),repmat(1,1,N));

% Inputs u(k), ..., u(k+N) (last one not used)
u = sdpvar(repmat(nu,1,N),repmat(1,1,N));

% Binary for PWA selection
d = binvar(2,1);

% Value functions
J = cell(1,N);

% Initialize value function at stage N
J{N} = 0;
t = sdpvar(nx+nu,1);
bounds(t,0,600);
for k = N-1:-1:1
    
    bounds(x{k},-5,5);
    bounds(u{k},-1,1);
    bounds(x{k+1},-5,5);
    
    % Feasible region
    F =     (-1 <= u{k}     <= 1);
    F = F + (-1 <= C*x{k}   <= 1);
    F = F + (-5 <= x{k}     <= 5);
    F = F + (-1 <= C*x{k+1} <= 1);
    F = F + (-5 <= x{k+1}   <= 5);

    % PWA Dynamics
    F = F + (implies(d(1),x{k+1} == (A*x{k}+B*u{k})));
    F = F + (implies(d(2),x{k+1} == (A*x{k}+pi*B*u{k})));
    F = F + (implies(d(1),x{k}(1) >= 0));
    F = F + (implies(d(2),x{k}(1) <= 0));
    F = F + (sum(d) == 1);
    
    F = F + (-t <= [x{k};u{k}] <= t) ;

   [mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,sum(t) + J{k+1},[],x{k},u{k});
 
end
mpsol{1} = mpt_removeOverlaps(mpsol{1})

% Compare
sysStruct.A{1} = A;
sysStruct.B{1} = B;
sysStruct.C{1} = C;
sysStruct.D{1} = [0];
sysStruct.A{2} = A;
sysStruct.B{2} = B*pi;
sysStruct.C{2} = C;
sysStruct.D{2} = [0];
sysStruct.guardX{1} = [-1 0];
sysStruct.guardU{1} = [0];
sysStruct.guardC{1} = [0];
sysStruct.guardX{2} = [1 0];
sysStruct.guardU{2} = [0];
sysStruct.guardC{2} = [0];

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
result =  fun1.join().compare(new.optimizer.join(), 'obj');
assert(result == 0)

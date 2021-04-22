function tests = test_mpt_pwa_3d_dp
tests = functiontests(localfunctions);

function test1(dummy)

% Prediction horizon
N = 4;

pwa3d
sysStruct.xmin = sysStruct.ymin;
sysStruct.xmax = sysStruct.ymax;
probStruct.R = 1;
probStruct.N=N-1;
probStruct.norm = 1;
probStruct.subopt_lev=0;
probStruct.P_N = zeros(3);

nx = 3;
nu = 1;

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
sysStruct.xmin = sysStruct.ymin;
sysStruct.xmax = sysStruct.ymax;
for k = N-1:-1:1
    % Feasible region
    t = sdpvar(nx+nu,1);
    bounds(x{k},sysStruct.xmin,sysStruct.xmax);
    bounds(u{k},sysStruct.umin,sysStruct.umax);
    bounds(x{k+1},sysStruct.xmin,sysStruct.ymax);
    bounds(t,0,600);

    F =     (sysStruct.umin <= u{k}     <= sysStruct.umax);
    F = F + (sysStruct.xmin <= x{k}     <= sysStruct.xmax);
    F = F + (sysStruct.xmin <= x{k+1}   <= sysStruct.xmax);
    F = F + (sysStruct.ymin <= sysStruct.C{1}*x{k}   <= sysStruct.ymax);
    F = F + (sysStruct.ymin <= sysStruct.C{1}*x{k+1} <= sysStruct.ymax);

    F = F + (-t <= [x{k};u{k}] <= t) ;

    % PWA Dynamics
    for i = 1:length(sysStruct.A)
        F = F + (implies(d(i),x{k+1} == sysStruct.A{i}*x{k}+sysStruct.B{i}*u{k}+sysStruct.f{i}));
        F = F + (implies(d(i),sysStruct.guardX{i}*x{k} <= sysStruct.guardC{i}));
    end
    F = F + (sum(d) == 1);

    % Compute value function for one step backwards
    [mpsol{k},sol{k},Uz{k},J{k}] = solvemp(F,sum(t) + J{k+1},[],x{k},u{k});
end
mpsol{1} = mpt_removeOverlaps(mpsol{1})


assert(length(mpsol{1}.Pn) == 182)

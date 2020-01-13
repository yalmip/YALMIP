function tests = test_gp_gpposybeam
tests = functiontests(localfunctions);

function test1(dummy)
N = 8;
w = sdpvar(N,1);
h = sdpvar(N,1);

% constants
wmin = .1; wmax = 100;
hmin = .1; hmax = 6;
Smin = 1/5; Smax = 5;
sigma_max = 1;
ymax = 10;
E = 1; F = 1;

% objective is the total volume of the beam
% obj = sum of (widths*heights*lengths) over each section
% (recall that the length of each segment is set to be 1)
obj = w'*h; 

% recursive formulation
v = sdpvar(N+1,1); y = sdpvar(N+1,1);
v(N+1,1) = 0; y(N+1,1) = 0;
for i = N:-1:1
  disp(['Processing recursion number: ' num2str(i)])
  v(i) = 12*(i-1/2)*F/(E*w(i)*h(i)^3) + v(i+1);
  y(i) = 6*(i-1/3)*F/(E*w(i)*h(i)^3)  + v(i+1) + y(i+1);
end

% constraint set
constr = [ ...
  wmin*ones(N,1) <= w; w <= wmax*ones(N,1);
  hmin*ones(N,1) <= h; h <= hmax*ones(N,1);
  Smin*ones(N,1) <= h./w; h./w <= Smax*ones(N,1);
  6*F*[1:N]'./(w.*(h.^2)) <= sigma_max*ones(N,1);
  y(1) <= ymax;
];

% solve GP and compute the optimal volume
sol  = optimize(constr,obj)

assert(sol.problem == 0);
assert(abs(value(obj)-42.39654132455499) <= 1e-3);

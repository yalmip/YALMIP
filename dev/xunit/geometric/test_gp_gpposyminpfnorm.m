function tests = test_gp_gpposyminpfnorm
tests = functiontests(localfunctions);

function test1(dummy)
sdpvar lambda b(4,1) s(3,1) v(4,1) c(2,1)

% constants
c_nom = [1 1]';
b_nom = [2 3 2 1]';
alpha = [1 1 1 1]'; beta  = [1 1 1 1]';
s_nom = [1 1 3]';
gamma = [1 1 1]'; delta = [1 1 1]';

% objective
obj = lambda;

% constraints
constr = [...
  % inequalities
  b'*v      <= lambda*v(1);
  s(1)*v(1) <= lambda*v(2);
  s(2)*v(2) <= lambda*v(3);
  s(3)*v(3) <= lambda*v(4);
  [0.5; 0.5] <= c; c <= [2; 2];
  % equalities
  b == b_nom.*((ones(4,1)*(c(1)/c_nom(1))).^alpha).*...
              ((ones(4,1)*(c(2)/c_nom(2))).^beta); 
  s == s_nom.*((ones(3,1)*(c(1)/c_nom(1))).^gamma).*...
              ((ones(3,1)*(c(2)/c_nom(2))).^delta);
];
constr=[constr, lambda >= 0, b >=0, s>=0, v>=0, c>=0];
% find the optimal eigenvalue
sol = optimize(constr,obj,sdpsettings('solver','mosek,gpposy,fmincon-geometric'));
obj
assert(sol.problem == 0);
assert(abs(value(obj) - 0.80406738656616) <= 1e-4);

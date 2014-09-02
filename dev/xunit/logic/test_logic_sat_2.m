function test_logic_sat_2

randn('seed',12345);
rand('seed',12345);

A1 = randn(8,2);
b1 = rand(8,1)*2-A1*[3;3];
A2 = randn(8,2);
b2 = rand(8,1)*2-A2*[-3;3];
A3 = randn(8,2);
b3 = rand(8,1)*2-A3*[3;-3];
A4 = randn(8,2);
b4 = rand(8,1)*2-A4*[-3;-3];

binvar inp1 inp2 inp3 inp4
F = true(inp1 | inp2 | inp3 | inp4);
x = sdpvar(2,1);
F = F + (iff(inp1,A1*x <= b1));
F = F + (iff(inp2,A2*x <= b2));
F = F + (iff(inp3,A3*x <= b3));
F = F + (iff(inp4,A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = solvesdp(F,-x(2));
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(x),[1.89084694511033   4.32311245145436]',1e-6);
mbg_asserttolequal(double([inp1 inp2 inp3 inp4]),[0 0 0 1]);

F = true(inp1 | inp2 | inp3 | inp4);
F = F + (inp1 == (A1*x <= b1));
F = F + (inp2 == (A2*x <= b2));
F = F + (inp3 == (A3*x <= b3));
F = F + (inp4 == (A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = solvesdp(F,-x(2));
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(x),[1.89084694511033   4.32311245145436]',1e-6);
mbg_asserttolequal(double([inp1 inp2 inp3 inp4]),[0 0 0 1]);

F = true(inp1 | inp2 | inp3 | inp4);
F = F + (implies(inp1,A1*x <= b1));
F = F + (implies(inp2,A2*x <= b2));
F = F + (implies(inp3,A3*x <= b3));
F = F + (implies(inp4,A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = solvesdp(F,-x(2));
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(x),[1.89084694511033   4.32311245145436]',1e-6);
mbg_asserttolequal(double([inp1 inp2 inp3 inp4]),[0 0 0 1]);

F = ( (A1*x <= b1) | (A2*x <= b2) | (A3*x <= b3) | (A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = solvesdp(F,-x(2));
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(x),[1.89084694511033   4.32311245145436]',1e-6);
mbg_asserttolequal(double([inp1 inp2 inp3 inp4]),[0 0 0 1]);

x = sdpvar(2,1);
bounds(x,-100,100);
F = ( (A1*x <= b1) | (A2*x <= b2) | (A3*x <= b3) | (A4*x <= b4));
F = F + (-100 <= x <= 100);
sol = solvesdp(F,-x(2));
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(x),[1.89084694511033   4.32311245145436]',1e-6);
mbg_asserttolequal(double([inp1 inp2 inp3 inp4]),[0 0 0 1]);

ii = sdpvar(1,1);
jj = sdpvar(1,1);
x = sdpvar(1,8);
p = [0 1 7 2 3 4 3 20];
solvesdp((-100 <= [x(:);ii;jj] <= 100) + (x == p)+(x([ii jj]) <= 3)+(ii~=jj),-ii-jj);
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(min(double([ii jj])),[5]);
mbg_asserttolequal(max(double([ii jj])),[7]);



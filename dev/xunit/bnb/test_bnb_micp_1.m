function tests = test_bnb_micp_1
tests = functiontests(localfunctions);

function test1(dummy)
randn('seed',123456);
n = 5;
%  65.63177408018771   6.95052852960134 -58.66230788287043 -70.86445052912004  58.37564024957875
%  17.58996579270265  65.53001366712473  18.43985990767680 -58.43579533652219 -72.32340882851518
% -46.45728349558073  15.51300078651673  68.32913865172714  18.82218749426777 -57.90310836270277
% -79.19194296790369 -53.92234161527504  21.23787385808123  65.53955212353939  17.48548000466040
%  51.25486176425121 -73.14184210768217 -50.35188380380193  19.83022901147247  67.22698658022694
P = toeplitz(randn(n,1)*100)+randn(n,n)*5;
Z = intvar(n,n,'toeplitz');
t = sdpvar(n,n,'full');
e = P(:)-Z(:);
ops = sdpsettings('solver','bnb','verbose',2);

F = (-t <= P-Z <= t);
obj = sum(sum(t));
sol = optimize(F,obj,ops);
assert(sol.problem == 0);
assert(abs(value(obj) - 66.18236738983525) <=  1e-4);

F = ([]);
obj = norm(e,1);
sol = optimize(F,obj,ops);
assert(sol.problem == 0);
assert(abs(value(obj) - 66.18236738983525) <= 1e-4);

obj = e'*e;
F = ([]);
sol = optimize(F,obj,ops);
assert(sol.problem == 0);
assert(abs(value(obj) - 3.352603490492911e+002) <= 1e-4);

t = sdpvar(1,1);
obj = t;
F = (cone(e,t));
sol = optimize(F,obj,ops);
assert(sol.problem == 0);
assert(abs(value(obj) - 18.31011603130778) <= 1e-4);

t = sdpvar(1,1);
obj = norm(e);
F = ([]);
sol = optimize(F,obj,ops);
assert(sol.problem == 0);
assert(abs(value(obj) - 18.31011603130778) <= 1e-4);

obj = t;
F = ([t e';e eye(length(e))]>=0);
sol = optimize(F,obj,ops);
assert(sol.problem == 0);
assert(abs(value(obj) - 3.352603420494530e+002) <= 2e-4);
function tests = test_fminconsdpvar
tests = functiontests(localfunctions);

function test1(dummy)
% Single LMI
yalmip('clear');
A = [1 0;0.4 1];
B = [0.4;0.08]; 
L = [1.9034 1.1501];  
Y = sdpvar(2,2);
F = [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0;
F = F+[(L*Y*L')<=1];
sol = optimize(F,-trace(Y),sdpsettings('solver','fmincon'));
assert(isequal(sol.problem,0));
assert(abs(value(trace(Y)-10.26039))<1e-3);

% Multiple LMI
F = [Y >=0, [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0];
F = F+[(L*Y*L')<=1];
sol = optimize(F,-trace(Y),sdpsettings('solver','fmincon'));
assert(isequal(sol.problem,0));
assert(abs(value(trace(Y)-10.26039))<1e-3);

% Some polynomial stuff should run at least
F = [Y'*Y >=0, [Y'*Y Y'*Y*(A-B*L)';(A-B*L)*Y'*Y Y'*Y]>=0];
F = F+[(L*Y'*Y*L')<=1];
sol = optimize(F,-trace(Y),sdpsettings('solver','fmincon'));
assert(isequal(sol.problem,0));
assert(abs(value(trace(Y)-3.5497))<1e-1);

% Some functional stuff should run at least
F = [Y >=0, [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0];
F = F+[exp((L*Y*L'))<=exp(1)];
sol = optimize(F,exp(-.1*trace(Y)),sdpsettings('solver','fmincon'));
assert(isequal(sol.problem,0));
assert(abs(value(trace(Y)-10.26039))<1e-1);

% Some functional stuff should run at least
F = [Y >=0, [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0];
F = F+[[exp(1)-exp((L*Y*L')) exp(1)-exp((L*Y*L'))]>=0];
sol = optimize(F,-trace(Y),sdpsettings('solver','fmincon','usex0',0));
assert(isequal(sol.problem,0));
assert(abs(value(trace(Y)-10.26039))<1e-1);

F = [Y >=0, [Y Y*(A-B*L)';(A-B*L)*Y Y]>=0];
F = F+[[exp(1)-exp((L*Y*L')) 0;0 exp(1)-exp((L*Y*L'))]>=0];
sol = optimize(F,-sqrtm(trace(Y)),sdpsettings('solver','fmincon','usex0',0));
assert(isequal(sol.problem,0));


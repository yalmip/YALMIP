function micp

yalmip('clear')
M = 3;
a = [3 4 5]';
N = 10;
x = intvar(1,2);gan(x,a(1:2));
x = [x intvar(1)];
assign(x,1);
sol = solvesdp(set(x >= 1)+set(sum(x)<=N),gan(x,a),sdpsettings('solver','bnb','usex0',1))
mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal( double(isequal(double(x), [1 2 2]) | isequal(double(x), [1 1 2])), 1);

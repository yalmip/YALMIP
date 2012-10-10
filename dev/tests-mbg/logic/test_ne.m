yalmip('clear')

binvar x y
solvesdp([0 <= [x,y] <= 1, x+y~=1],(x+y-1)^2-10*x)
mbg_asserttolequal(double(x),1, 1e-4);
mbg_asserttolequal(double(y),1, 1e-4);

binvar x y
solvesdp([0 <= [x,y] <= 1, x+y~=1],(x+y-1)^2)
mbg_asserttolequal(double((x+y-1)^2),1, 1e-4);
function gp

sdpvar x y
css=set(x>=1)+set(y>=1)+set(x/y<=4)+set(y<=8);
css=css+set(x^2/y==1.5);
obj = x+y/x;
sol = solvesdp(css,obj)

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(obj),2.04124145231932,1e-4);
mbg_asserttolequal(double([x y]), [ 1.22474487139159   1.00000000000000], 1e-4);

sdpvar x y z
css=set(x>=1)+set(y>=1)+set(x/y<=4)+set(y<=8) + set(x*z == 10) %+ set(1<z<16);
css=css+set(x^2/y==1.5);
obj = x+y/x;
sol = solvesdp(css,obj)

mbg_asserttolequal(sol.problem,0);
mbg_asserttolequal(double(obj),2.04124145231932,1e-4);
mbg_asserttolequal(double([x y z]), [ 1.22474487139159   1.00000000000000 8.16496580927726], 1e-4);


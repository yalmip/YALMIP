function test_sdpvar_norm

% To improve performance, we don't introduce auxilliary variables for
% elements which are fixed
sdpvar x y
assign(x,1);
mbg_asserttolequal(double(norm([x;-3],1)),4);
solvesdp(norm([x;-3],1)<=y,y);
mbg_asserttolequal(double(y),3,1e-4)
mbg_asserttolequal(double(x),0,1e-4)

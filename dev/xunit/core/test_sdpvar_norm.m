function test_sdpvar_norm

% To improve performance, we don't introduce auxilliary variables for
% elements which are fixed
sdpvar x y
assign(x,1);
assertElementsAlmostEqual(double(norm([x;-3],1)),4,'absolute',1e-8);
solvesdp(norm([x;-3],1)<=y,y);
assertElementsAlmostEqual(double(y),3,'absolute',1e-8)
assertElementsAlmostEqual(double(x),0,'absolute',1e-8)

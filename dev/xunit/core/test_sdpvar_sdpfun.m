function test_sdpvar_sdpfun

yalmip('clear');
sdpvar x
assign(x,0.1);
ops = sdpsettings('debug',1,'fmincon.algorithm','sqp','usex0',1);
sol = solvesdp(sdpfun(x,1,'mytestOLD') >= 0,x,ops)
assertTrue(sol.problem == 0);

yalmip('clear');
sdpvar x
assign(x,0.1);
ops = sdpsettings('debug',1,'fmincon.algorithm','sqp','usex0',1);
sol = solvesdp(sdpfun(1,x,'mytestNEW') >= 0,x,ops)
assertTrue(sol.problem == 0);


function tests = test_sdpvar_sdpfun
tests = functiontests(localfunctions);

function test1(dummy)
yalmip('clear');
sdpvar x
assign(x,0.1);
ops = sdpsettings('debug',1,'fmincon.algorithm','sqp','usex0',1);
sol = optimize(sdpfun(x,1,'mytestOLD') >= 0,x,ops)
assert(sol.problem == 0);

yalmip('clear');
sdpvar x
assign(x,0.1);
ops = sdpsettings('debug',1,'fmincon.algorithm','sqp','usex0',1);
sol = optimize(sdpfun(1,x,'mytestNEW') >= 0,x,ops)
assert(sol.problem == 0);


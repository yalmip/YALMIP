function sdpvar_toeplitz

c = sdpvar(4,1);
assign(c,[2 4 6 8]');
ok = all(all(double(toeplitz(c,c))-toeplitz(double(c),double(c)) == 0));
mbg_asserttrue(ok)

ok = all(all(double(toeplitz(c))-toeplitz(double(c)) == 0));
mbg_asserttrue(ok)

r = sdpvar(4,1);
assign(r,-[2 4 6 8]');
ok = all(all(double(toeplitz(c,r))-toeplitz(double(c),double(r)) == 0));
mbg_asserttrue(ok)

c = sdpvar(4,1);
assign(c,[2 4 6 8]');
r = sdpvar(2,1);
assign(r,-[2 4]');
ok = all(all(double(toeplitz(c,r))-toeplitz(double(c),double(r)) == 0));
mbg_asserttrue(ok)

ok = all(all(double(toeplitz(r,c'))-toeplitz(double(r),double(c')) == 0));
mbg_asserttrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]');
ok = all(all(double(toeplitz(c,r))-toeplitz(double(c),double(r))==0));
mbg_asserttrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]');
ok = all(all(double(toeplitz(r,c))-toeplitz(double(r),double(c))==0));
mbg_asserttrue(ok)

c = sdpvar(4,1,'full','complex');
assign(c,[2 4 6 8]'+sqrt(-1)*[5 4 3 2]');
ok = all(all(double(toeplitz(c,c))-toeplitz(double(c),double(c)) == 0));
mbg_asserttrue(ok)

ok = all(all(double(toeplitz(c))-toeplitz(double(c)) == 0));
mbg_asserttrue(ok)

c = sdpvar(4,1,'full','complex');
assign(c,[2 4 6 8]'+sqrt(-1)*[5 4 3 2]');
r = sdpvar(4,1);
assign(r,-[2 4 6 8]');
ok = all(all(double(toeplitz(c,r))-toeplitz(double(c),double(r)) == 0));
mbg_asserttrue(ok)

r = sdpvar(2,1);
assign(r,-[2 4]'+sqrt(-1));
ok = all(all(double(toeplitz(c,r))-toeplitz(double(c),double(r)) == 0));
mbg_asserttrue(ok)

ok = all(all(double(toeplitz(r,c'))-toeplitz(double(r),double(c')) == 0));
mbg_asserttrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]'*sqrt(-1));
ok = all(all(double(toeplitz(c,r))-toeplitz(double(c),double(r))==0));
mbg_asserttrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]*sqrt(-1));
assign(r,[1 2 3 4]');
ok = all(all(double(toeplitz(r,c))-toeplitz(double(r),double(c))==0));
mbg_asserttrue(ok)




function sdpvar_hankel

c = sdpvar(4,1);
assign(c,[2 4 6 8]');
ok = all(all(double(hankel(c,c))-hankel(double(c),double(c)) == 0));
assertTrue(ok)

ok = all(all(double(hankel(c))-hankel(double(c)) == 0));
assertTrue(ok)

r = sdpvar(4,1);
assign(r,-[2 4 6 8]');
ok = all(all(double(hankel(c,r))-hankel(double(c),double(r)) == 0));
assertTrue(ok)

c = sdpvar(4,1);
assign(c,[2 4 6 8]');
r = sdpvar(2,1);
assign(r,-[2 4]');
ok = all(all(double(hankel(c,r))-hankel(double(c),double(r)) == 0));
assertTrue(ok)

ok = all(all(double(hankel(r,c'))-hankel(double(r),double(c')) == 0));
assertTrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]');
ok = all(all(double(hankel(c,r))-hankel(double(c),double(r))==0));
assertTrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]');
ok = all(all(double(hankel(r,c))-hankel(double(r),double(c))==0));
assertTrue(ok)

c = sdpvar(4,1,'full','complex');
assign(c,[2 4 6 8]'+sqrt(-1)*[5 4 3 2]');
ok = all(all(double(hankel(c,c))-hankel(double(c),double(c)) == 0));
assertTrue(ok)

ok = all(all(double(hankel(c))-hankel(double(c)) == 0));
assertTrue(ok)

c = sdpvar(4,1,'full','complex');
assign(c,[2 4 6 8]'+sqrt(-1)*[5 4 3 2]');
r = sdpvar(4,1);
assign(r,-[2 4 6 8]');
ok = all(all(double(hankel(c,r))-hankel(double(c),double(r)) == 0));
assertTrue(ok)

r = sdpvar(2,1);
assign(r,-[2 4]'+sqrt(-1));
ok = all(all(double(hankel(c,r))-hankel(double(c),double(r)) == 0));
assertTrue(ok)

ok = all(all(double(hankel(r,c'))-hankel(double(r),double(c')) == 0));
assertTrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]'*sqrt(-1));
ok = all(all(double(hankel(c,r))-hankel(double(c),double(r))==0));
assertTrue(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]*sqrt(-1));
assign(r,[1 2 3 4]');
ok = all(all(double(hankel(r,c))-hankel(double(r),double(c))==0));
assertTrue(ok)













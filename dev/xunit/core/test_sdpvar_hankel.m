function tests = test_sdpvar_hankel
tests = functiontests(localfunctions);

function test1(dummy)

c = sdpvar(4,1);
assign(c,[2 4 6 8]');
ok = all(all(value(hankel(c,c))-hankel(value(c),value(c)) == 0));
assert(ok)

ok = all(all(value(hankel(c))-hankel(value(c)) == 0));
assert(ok)

r = sdpvar(4,1);
assign(r,-[2 4 6 8]');
ok = all(all(value(hankel(c,r))-hankel(value(c),value(r)) == 0));
assert(ok)

c = sdpvar(4,1);
assign(c,[2 4 6 8]');
r = sdpvar(2,1);
assign(r,-[2 4]');
ok = all(all(value(hankel(c,r))-hankel(value(c),value(r)) == 0));
assert(ok)

ok = all(all(value(hankel(r,c'))-hankel(value(r),value(c')) == 0));
assert(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]');
ok = all(all(value(hankel(c,r))-hankel(value(c),value(r))==0));
assert(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]');
ok = all(all(value(hankel(r,c))-hankel(value(r),value(c))==0));
assert(ok)

c = sdpvar(4,1,'full','complex');
assign(c,[2 4 6 8]'+sqrt(-1)*[5 4 3 2]');
ok = all(all(value(hankel(c,c))-hankel(value(c),value(c)) == 0));
assert(ok)

ok = all(all(value(hankel(c))-hankel(value(c)) == 0));
assert(ok)

c = sdpvar(4,1,'full','complex');
assign(c,[2 4 6 8]'+sqrt(-1)*[5 4 3 2]');
r = sdpvar(4,1);
assign(r,-[2 4 6 8]');
ok = all(all(value(hankel(c,r))-hankel(value(c),value(r)) == 0));
assert(ok)

r = sdpvar(2,1);
assign(r,-[2 4]'+sqrt(-1));
ok = all(all(value(hankel(c,r))-hankel(value(c),value(r)) == 0));
assert(ok)

ok = all(all(value(hankel(r,c'))-hankel(value(r),value(c')) == 0));
assert(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]);
assign(r,[1 2 3 4]'*sqrt(-1));
ok = all(all(value(hankel(c,r))-hankel(value(c),value(r))==0));
assert(ok)

c = sdpvar(2,3);
r = sdpvar(4,1);
assign(c,[1 2 3;4 5 6]*sqrt(-1));
assign(r,[1 2 3 4]');
ok = all(all(value(hankel(r,c))-hankel(value(r),value(c))==0));
assert(ok)
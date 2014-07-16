function test_rcone

x = sdpvar(1);
y = sdpvar(1);
z = sdpvar(3,1);

solvesdp([z>=1, rcone(z,x,y)],x+y);
assertTrue(abs(double(x*y)-1.5) <= 1e-4);
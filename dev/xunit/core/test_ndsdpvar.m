function test_ndsdpvar

x = sdpvar(1,2,3);
p = randn(1,2,3);
assign(x,p)

assertTrue(isequal(double(x),p))
assertTrue(isequal(double(x(:)),p(:)))

assertTrue(isequal(double(sum(x)),sum(p)))
assertTrue(isequal(double(sum(x,1)),sum(p,1)))
assertTrue(isequal(double(sum(x,2)),sum(p,2)))
assertTrue(isequal(double(sum(x,3)),sum(p,3)))

assertTrue(isequal(double(diff(x)),diff(p)))
assertTrue(isequal(double(diff(x,1)),diff(p,1)))
assertVectorsAlmostEqual(double(diff(x,2)),diff(p,2))
assertVectorsAlmostEqual(double(diff(x,3)),diff(p,3))

assertTrue(isequal(double(diff(x)),diff(p)))
assertTrue(isequal(double(diff(x,1,1)),diff(p,1,1)))
assertTrue(isequal(double(diff(x,2,1)),diff(p,2,1)))
assertTrue(isequal(double(diff(x,3,1)),diff(p,3,1)))

assertTrue(isequal(double(diff(x)),diff(p)))
assertTrue(isequal(double(diff(x,1,2)),diff(p,1,2)))
assertTrue(isequal(double(diff(x,2,2)),diff(p,2,2)))
assertTrue(isequal(double(diff(x,3,2)),diff(p,3,2)))

x = sdpvar(1,2,3);
y = sdpvar(1,2,3);
p1 = randn(1,2,3);
p2 = randn(1,2,3);
assign(x,p1)
assign(y,p2)
z =  x - y;
assertTrue(isequal(double(z(:)),p1(:)-p2(:)))

x = sdpvar(1,2,3);
y = sdpvar(1,1);
p1 = randn(1,2,3);
p2 = randn(1);
assign(x,p1)
assign(y,p2)
z =  x - y;
assertTrue(isequal(double(z(:)),p1(:)-p2(:)))

X = randn(2,3,4);
Xv = sdpvar(2,3,4);
assign(Xv,X);
assertTrue(norm(double(max(Xv,[],3))-max(X,[],3)) <= 1e-6);
s1 = double(max(Xv,[],2));
s2 = max(X,[],2);
assertTrue(norm(s1(:)-s2(:)) <= 1e-6);
s1 = double(max(Xv,[],1));
s2 = max(X,[],1);
assertTrue(norm(s1(:)-s2(:)) <= 1e-6);
X = randn(4,2,3,2);
Xv = sdpvar(4,2,3,2);
assign(Xv,X);
s1 = double(max(Xv,[],4));
s2 = max(X,[],4);
assertTrue(norm(s1(:)-s2(:)) <= 1e-6);
s1 = double(max(Xv,[],3));
s2 = max(X,[],3);
assertTrue(norm(s1(:)-s2(:)) <= 1e-6);
s1 = double(max(Xv,[],2));
s2 = max(X,[],2);
assertTrue(norm(s1(:)-s2(:)) <= 1e-6);
s1 = double(max(Xv,[],1));
s2 = max(X,[],1);
assertTrue(norm(s1(:)-s2(:)) <= 1e-6);




function r = isequal(a,b)
r = (norm(a(:)-b(:)) < 1e-12) & all(size(a) == size(b));
function tests = test_ndsdpvar
tests = functiontests(localfunctions);

function test1(dummy)

x = sdpvar(1,2,3);
p = randn(1,2,3);
assign(x,p)

assert(isequal(value(x),p))
assert(isequal(value(x(:)),p(:)))

assert(isequal(value(sum(x)),sum(p)))
assert(isequal(value(sum(x,1)),sum(p,1)))
assert(isequal(value(sum(x,2)),sum(p,2)))
assert(isequal(value(sum(x,3)),sum(p,3)))

assert(isequal(value(diff(x)),diff(p)))
assert(isequal(value(diff(x,1)),diff(p,1)))
assert(norm(sum(value(diff(x,2))-diff(p,2))) <= 1e-4)
assert(norm(sum(value(diff(x,3))-diff(p,3))) <= 1e-4)

assert(isequal(value(diff(x)),diff(p)))
assert(isequal(value(diff(x,1,1)),diff(p,1,1)))
assert(isequal(value(diff(x,2,1)),diff(p,2,1)))
assert(isequal(value(diff(x,3,1)),diff(p,3,1)))

assert(isequal(value(diff(x)),diff(p)))
assert(isequal(value(diff(x,1,2)),diff(p,1,2)))
assert(isequal(value(diff(x,2,2)),diff(p,2,2)))
assert(isequal(value(diff(x,3,2)),diff(p,3,2)))

x = sdpvar(1,2,3);
y = sdpvar(1,2,3);
p1 = randn(1,2,3);
p2 = randn(1,2,3);
assign(x,p1)
assign(y,p2)
z =  x - y;
assert(isequal(value(z(:)),p1(:)-p2(:)))

x = sdpvar(1,2,3);
y = sdpvar(1,1);
p1 = randn(1,2,3);
p2 = randn(1);
assign(x,p1)
assign(y,p2)
z =  x - y;
assert(isequal(value(z(:)),p1(:)-p2(:)))

X = randn(2,3,4);
Xv = sdpvar(2,3,4);
assign(Xv,X);
assert(norm(value(max(Xv,[],3))-max(X,[],3)) <= 1e-6);
s1 = value(max(Xv,[],2));
s2 = max(X,[],2);
assert(norm(s1(:)-s2(:)) <= 1e-6);
s1 = value(max(Xv,[],1));
s2 = max(X,[],1);
assert(norm(s1(:)-s2(:)) <= 1e-6);
X = randn(4,2,3,2);
Xv = sdpvar(4,2,3,2);
assign(Xv,X);
s1 = value(max(Xv,[],4));
s2 = max(X,[],4);
assert(norm(s1(:)-s2(:)) <= 1e-6);
s1 = value(max(Xv,[],3));
s2 = max(X,[],3);
assert(norm(s1(:)-s2(:)) <= 1e-6);
s1 = value(max(Xv,[],2));
s2 = max(X,[],2);
assert(norm(s1(:)-s2(:)) <= 1e-6);
s1 = value(max(Xv,[],1));
s2 = max(X,[],1);
assert(norm(s1(:)-s2(:)) <= 1e-6);

function r = isequal(a,b)
r = (norm(a(:)-b(:)) < 1e-12) & all(size(a) == size(b));
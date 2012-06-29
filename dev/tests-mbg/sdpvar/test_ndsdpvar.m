function test_ndsdpvar

x = sdpvar(1,2,3);
p = randn(1,2,3);
assign(x,p)

mbg_asserttrue(isequal(double(x),p))
mbg_asserttrue(isequal(double(x(:)),p(:)))

mbg_asserttrue(isequal(double(sum(x)),sum(p)))
mbg_asserttrue(isequal(double(sum(x,1)),sum(p,1)))
mbg_asserttrue(isequal(double(sum(x,2)),sum(p,2)))
mbg_asserttrue(isequal(double(sum(x,3)),sum(p,3)))

mbg_asserttrue(isequal(double(diff(x)),diff(p)))
mbg_asserttrue(isequal(double(diff(x,1)),diff(p,1)))
mbg_asserttrue(isequal(double(diff(x,2)),diff(p,2)))
mbg_asserttrue(isequal(double(diff(x,3)),diff(p,3)))

mbg_asserttrue(isequal(double(diff(x)),diff(p)))
mbg_asserttrue(isequal(double(diff(x,1,1)),diff(p,1,1)))
mbg_asserttrue(isequal(double(diff(x,2,1)),diff(p,2,1)))
mbg_asserttrue(isequal(double(diff(x,3,1)),diff(p,3,1)))

mbg_asserttrue(isequal(double(diff(x)),diff(p)))
mbg_asserttrue(isequal(double(diff(x,1,2)),diff(p,1,2)))
mbg_asserttrue(isequal(double(diff(x,2,2)),diff(p,2,2)))
mbg_asserttrue(isequal(double(diff(x,3,2)),diff(p,3,2)))

x = sdpvar(1,2,3);
y = sdpvar(1,2,3);
p1 = randn(1,2,3);
p2 = randn(1,2,3);
assign(x,p1)
assign(y,p2)
z =  x - y;
mbg_asserttrue(isequal(double(z(:)),p1(:)-p2(:)))

x = sdpvar(1,2,3);
y = sdpvar(1,1);
p1 = randn(1,2,3);
p2 = randn(1);
assign(x,p1)
assign(y,p2)
z =  x - y;
mbg_asserttrue(isequal(double(z(:)),p1(:)-p2(:)))



function r = isequal(a,b)
r = (norm(a(:)-b(:)) < 1e-12) & all(size(a) == size(b));
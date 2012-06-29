function test_diff


% border cases
mbg_asserttrue(isempty(diff(sdpvar(1,1))));

U = sdpvar(1,1);
mbg_asserttrue(isequal(0,diff([U U])));
mbg_asserttrue(isequal(0,diff([U;U])));

mbg_asserttrue(isequal(U - diff([0 U]),0));
mbg_asserttrue(isequal(-U - diff([U 0]),0));

U = sdpvar(2,4);
r = randn(2,4);
assign(U,r)
mbg_asserttrue(norm(double(diff(U)-diff(r)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1)-diff(r,1)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2)-diff(r,2)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1,2)-diff(r,1,2)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2,2)-diff(r,2,2)))<1e-12)

U = sdpvar(4,2);
r = randn(4,2);
assign(U,r)
mbg_asserttrue(norm(double(diff(U)-diff(r)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1)-diff(r,1)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2)-diff(r,2)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1,2)-diff(r,1,2)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2,2)-diff(r,2,2)))<1e-12)

U = sdpvar(4,4);
r = randn(4,4);r=r+r';
assign(U,r)
mbg_asserttrue(norm(double(diff(U)-diff(r)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1)-diff(r,1)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2)-diff(r,2)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1,2)-diff(r,1,2)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2,2)-diff(r,2,2)))<1e-12)

U = sdpvar(4,1);
r = randn(4,1);
assign(U,r)
mbg_asserttrue(norm(double(diff(U)-diff(r)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1)-diff(r,1)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2)-diff(r,2)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1,2)-diff(r,1,2)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2,2)-diff(r,2,2)))<1e-12)

U = sdpvar(4,3,'fu','co');
r = randn(4,3)+randn(4,3)*sqrt(-1);
assign(U,r)
mbg_asserttrue(norm(double(diff(U)-diff(r)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1)-diff(r,1)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2)-diff(r,2)))<1e-12)

mbg_asserttrue(norm(double(diff(U,1,2)-diff(r,1,2)))<1e-12)
mbg_asserttrue(norm(double(diff(U,2,2)-diff(r,2,2)))<1e-12)


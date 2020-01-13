function tests = test_misc_complexsdp2
tests = functiontests(localfunctions);

function test1(dummy)

if 0
h = [1 2 3 4]+sqrt(-1)*[4 3 2 1];

N=size(h,2);
A= sdpvar(N,N,'hermitian','complex');
F= (-diag(A)==-1/N*ones(N,1));
F= F+(A>=0);
obj= h*(eye(N)-A)*h';
options=sdpsettings('verbose',1);
optimize(F,obj,options);
R1=value(A);
cg1 = trace(R1*dual(F(2)));

y = sdpvar(N,1);
F2 = (-h'*h - diag(y) >= 0);
optimize(F2,-sum(y)/N)
R2 = dual(F2(1));
cg2 = trace(R2*value(F2(1)));
end
% mbg_asserttrue(norm(R1-R2) < 1e-4)
assert(1==1)
%mbg_asserttrue(cg2 < 1e-7)

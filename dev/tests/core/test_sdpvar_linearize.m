function tests = test_linearize
tests = functiontests(localfunctions);

function test1(testCase) % univariate_scalar_poly_lin
x = sdpvar();
x0 = 5;
assign(x,x0);
p = x^5+4*x^4+3*x^3+2*x^2+x-1;
dp0 = value(5*x^4+16*x^3+9*x^2+4*x+1);
p0 = value(p);
pLin = linearize(p);
assign(x,6)
testCase.assertTrue(norm(value(pLin)-value(p0+dp0'*(x-x0)))<1e-10);

function test2(testCase) % multivariate_scalar_poly_lin
n = 10;
x = sdpvar(n,1);
p = ((x-ones(n,1))'*(x-ones(n,1)))^2;
dp = 4*(x-ones(n,1))'*(x-ones(n,1))*(x-ones(n,1));
x0 = (1:n)';
assign(x,x0);
p0 = value(p);
dp0 = value(dp);
pLin = linearize(p);
assign(x,(1:n)*2);
testCase.assertTrue(norm(value(pLin)-value(p0+dp0'*(x-x0)))<1e-10);


function test3(testCase) % multivariate_matrix_poly_lin
n = 10;
m = 5;
X = sdpvar(n,n,'full');
B = randn(n,m);
C = randn(n,m);
p = B'*X'*X*C;
h = @(x) x+x';
X0 = randn(n,n);
assign(X,X0);
dp = cell(m);
for i = 1:m
	for j = 1:m
		dp{i,j} = value(X*h(B(:,i)*C(:,j)'));
	end
end
p0 = value(p);
pLin = linearize(p);

assign(X,randn(n,n));
dX = value(X-X0);
dpdX = zeros(m,m);
for i = 1:m
	for j = 1:m
		dpdX(i,j) = trace(dp{i,j}'*dX);
	end
end
pLinCorrect = p0 + dpdX;

testCase.assertTrue(norm(pLinCorrect-value(pLin))<=1e-10);
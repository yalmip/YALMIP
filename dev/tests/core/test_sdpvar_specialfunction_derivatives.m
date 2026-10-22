function tests = test_sdpvar_specialfunction_derivatives
tests = functiontests(localfunctions);

function test_asech(testCase)
checkDerivative(testCase,@asech,[0.1 0.3 0.5 0.9]);

function test_asec(testCase)
checkDerivative(testCase,@asec,[-3 -1.5 1.5 3]);

function test_acsc(testCase)
checkDerivative(testCase,@acsc,[-3 -1.5 1.5 3]);

function test_coth(testCase)
checkDerivative(testCase,@coth,[-3 -0.3 0.3 3]);

function test_atanh(testCase)
checkDerivative(testCase,@atanh,[-0.9 -0.5 0 0.5 0.9]);

function checkDerivative(testCase,f,points)
x = sdpvar(1);
y = f(x);
e = yalmip('extstruct',getvariables(y));
[~,operator] = feval(e.fcn,'graph',e.var,e.arg{1:end-1});

% Compare the registered callback with the numeric function, including
% both real branches where applicable and scalar/row/column inputs.
inputs = {points(1),points,points'};
for k = 1:numel(inputs)
    z = inputs{k};
    h = 1e-6*max(1,abs(z));
    expected = (f(z+h)-f(z-h))./(2*h);
    actual = operator.derivative(z);
    testCase.assertTrue(isreal(actual));
    testCase.assertTrue(all(isfinite(actual(:))));
    testCase.verifyEqual(actual,expected,'RelTol',1e-5,'AbsTol',1e-8);
end

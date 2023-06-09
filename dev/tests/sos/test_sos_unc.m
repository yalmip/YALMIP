function tests = test_sos_unc
tests = functiontests(localfunctions);

function test1(testCase)

sdpvar x y a t 
gamma = sdpvar(1);

p = a*x^4 + y^4+x*y + 1+gamma;
sol = solvesos((uncertain(gamma))+(-1/2 <= gamma <= 1/2)+(sos(p - t))+(ismember(a,[3 4 5]))+(4.6 >= abs(a) >= 3.5),-t,sdpsettings('verbose',0),[a t])

testCase.assertTrue(sol.problem == 0)
testCase.assertTrue(abs(value(t) - 0.4375) <= 1e-5)
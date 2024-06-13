function tests = test_operator_optimizer11
tests = functiontests(localfunctions);

function test1(testCase)

% Tests a regression bug that made expandmodel flawed. Basically, when
% optmizer generates the model, it constraints the parametric variables to
% pi. However, these constraints can not be used to tighten the model since
% they are completely artificial.

sdpvar x y
P = optimizer([x <= sin(y)*cos(y)],x^2,sdpsettings('solver','+quadprog'),y,x)
testCase.assertTrue(abs(P{-pi/3}--.433) <= 1e-3)

function test2(testCase)
sdpvar y x
P = optimizer([x <= sin(y)*cos(y)],x^2,sdpsettings('solver','+quadprog'),y,x)
testCase.assertTrue(abs(P{-pi/3}--.433) <= 1e-3)

function test3(testCase)
sdpvar y x z
P = optimizer([x <= sin(y)*cos(z)],x^2,sdpsettings('solver','+quadprog'),[y z],x)
testCase.assertTrue(abs(P{[-pi/3 -pi/5]}--.7006) <= 1e-3)

function test4(testCase)
sdpvar y x z
P = optimizer([x <= sin(y)*cos(z)],x^2,sdpsettings('solver','+quadprog'),[z y],x)
testCase.assertTrue(abs(P{[-pi/5 -pi/3]}--.7006) <= 1e-3)

function test5(testCase)
sdpvar x y  z
P = optimizer([x <= sin(y)*cos(z)],x^2,sdpsettings('solver','+quadprog'),[z y],x)
testCase.assertTrue(abs(P{[-pi/5 -pi/3]}--.7006) <= 1e-3)


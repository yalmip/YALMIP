function tests = test_bnb_dual_certificate
tests = functiontests(localfunctions);
end

function test_psd_factor_gives_original_model_lower_bound(testCase)
% Minimize x subject to 0<=x<=1 and diag(x,1)>=0. The exact bound is 0.
At = sparse([-1;1;-1;0;0;0]);
c = [0;1;0;0;0;1];
b = -1;
K.l = 2; K.s = 2;
B = [1 0;0 0];
[bound,details] = bnbCertifySDPDualLPBound(At,b,c,K,B,[0;0]);
testCase.verifyEqual(details.reason,'finite');
testCase.verifyGreaterThan(bound,-1e-9);
testCase.verifyLessThanOrEqual(bound,0);
end

function test_unbounded_reduced_cost_fails_closed(testCase)
% The Lagrangian -x is unbounded below when only x>=0 remains.
At = sparse([-1;-1;0;0;0]);
c = [0;0;0;0;1];
b = 0;
K.l = 1; K.s = 2;
B = [1 0;0 0];
[bound,details] = bnbCertifySDPDualLPBound(At,b,c,K,B,0);
testCase.verifyEqual(bound,-inf);
testCase.verifyEqual(details.reason,'unbounded residual direction');
end

function test_explicit_box_and_perturbed_dual(testCase)
% Arbitrary nonnegative LP multiplier and PSD factor must stay below every
% feasible objective, even without stationarity.
At = sparse([-1;1;-1;0;0;0]);
c = [0;1;0;0;0;1];
b = -1;
K.l = 2; K.s = 2;
B = [0.7 0.2;0.1 0.3];
y = [0.4;0.2];
[bound,details] = bnbCertifySDPDualLPBound(At,b,c,K,B,y, ...
    sparse(0,1),zeros(0,1),0,1);
testCase.verifyEqual(details.reason,'finite');
testCase.verifyLessThanOrEqual(bound,0);
for x = linspace(0,1,21)
    testCase.verifyLessThanOrEqual(bound,x+1e-12);
end
end

function test_random_psd_factors_never_raise_the_known_optimum(testCase)
rng(17,'twister');
K.l = 2; K.s = 2;
for trial = 1:50
    D = randn(2);
    D = 0.15*(D+D')/2;
    At = sparse([-1;1;-D(:)]);
    c = [0;1;1;0;0;1];
    objectiveCoefficient = randn;
    b = -objectiveCoefficient;
    B = randn(2);
    y = abs(randn(2,1));
    [bound,details] = bnbCertifySDPDualLPBound(At,b,c,K,B,y);
    testCase.verifyEqual(details.reason,'finite');
    testCase.verifyLessThanOrEqual(bound,min(0,objectiveCoefficient)+1e-10);
end
end

function test_rectangular_factor_and_extreme_scale_fallback(testCase)
At = sparse([-1;1;-1;0;0;0]);
c = [0;1;0;0;0;1];
K.l = 2; K.s = 2;
[bound,details] = bnbCertifySDPDualLPBound(At,-1,c,K,[0.5 0.1],[0;0]);
testCase.verifyEqual(details.reason,'finite');
testCase.verifyLessThanOrEqual(bound,0);
[bound,details] = bnbCertifySDPDualLPBound(At,-1,c,K,[1e-160 0],[0;0]);
testCase.verifyEqual(bound,-inf);
testCase.verifyEqual(details.reason,'extreme input scale');
end

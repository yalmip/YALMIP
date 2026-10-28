function tests = test_bnb_objective_spectrum
tests = functiontests(localfunctions);
end

function test_two_cost_classes_strengthen_a_fractional_bound(testCase)
p = model([0.5;0.5;sqrt(2)],zeros(3,1),ones(3,1),0);
[s,info] = bnbInitObjectiveSpectrum(p,1+sqrt(2));
testCase.verifyEqual(info.reason,'active');
bound = bnbRoundObjectiveSpectrum(s,1.1);
testCase.verifyGreaterThan(bound,1.4);
testCase.verifyLessThanOrEqual(bound,sqrt(2));
end

function test_exhaustive_signed_integer_domains_and_near_equal_costs(testCase)
rng(1729,'twister');
digits = dec2base(0:3^6-1,3,6)-'0';
for trial = 1:40
    atoms = [0.25;sqrt(0.5)];
    c = atoms(randi(2,6,1)).*sign(rand(6,1)-0.5);
    if mod(trial,2), c = c+1e-12*randn(6,1); end
    lo = randi([-3 2],6,1);
    p = model(c,lo,lo+2,randn);
    points = digits+lo';
    objectives = sort(p.f+points*c);
    upper = objectives(randi([400 numel(objectives)]));
    s = bnbInitObjectiveSpectrum(p,upper);
    testCase.assertNotEmpty(s);
    for j = 1:20
        lower = objectives(1)+(upper-objectives(1))*rand;
        oracle = min(objectives(objectives>=lower));
        bound = bnbRoundObjectiveSpectrum(s,lower);
        testCase.verifyGreaterThanOrEqual(bound,lower);
        testCase.verifyLessThanOrEqual(bound,oracle+1e-11);
    end
end
end

function test_continuous_nonlinear_and_large_spectra_are_skipped(testCase)
p = model([0.5;sqrt(2)],zeros(2,1),ones(2,1),0);
p.integral_variables = 1;
testCase.verifyEmpty(bnbInitObjectiveSpectrum(p,3));
p.integral_variables = 1:2; p.Q(1,1) = 1;
testCase.verifyEmpty(bnbInitObjectiveSpectrum(p,3));
p.Q(:) = 0; p.IntegerCost = true;
testCase.verifyEmpty(bnbInitObjectiveSpectrum(p,3));
p = model(sqrt([2;3;5;7]),zeros(4,1),100*ones(4,1),0);
testCase.verifyEmpty(bnbInitObjectiveSpectrum(p,1000));
testCase.verifyEqual(bnbRoundObjectiveSpectrum([],2.5),2.5);
end

function test_fixed_objective_terms_and_zero_cost_continuous_variables(testCase)
p = model([-sqrt(2);0.5;0],[2;0;-inf],[2;3;inf],-0.3);
p.integral_variables = 1:2;
s = bnbInitObjectiveSpectrum(p,-0.3-2*sqrt(2)+1.5);
testCase.assertNotEmpty(s);
lower = -0.3-2*sqrt(2)+0.7;
bound = bnbRoundObjectiveSpectrum(s,lower);
testCase.verifyGreaterThan(bound,lower+0.2);
testCase.verifyLessThanOrEqual(bound,-0.3-2*sqrt(2)+1);
end

function p = model(c,lo,hi,f)
n = numel(c);
p = struct('c',c,'f',f,'Q',sparse(n,n),'evalMap',[], ...
    'variabletype',zeros(n,1),'K',struct('m',0),'IntegerCost',false, ...
    'integral_variables',1:n,'lb',lo,'ub',hi);
end

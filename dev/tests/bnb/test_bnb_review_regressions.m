function tests = test_bnb_review_regressions
tests = functiontests(localfunctions);

function test_binary_product_small_tables(testCase)
for n = 1:3
    p.binaryProduct = reshape(1:3*n,3,n)';
    p.binaryProduct = p.binaryProduct(:,[3 1 2]);
    p.feasible = 1;
    for fixed = 1:4
        p.lb = zeros(3*n,1); p.ub = ones(3*n,1);
        if fixed == 1
            p.lb(1:2) = 1;
        elseif fixed == 2
            p.ub(1) = 0;
        elseif fixed == 3
            p.lb(1) = 1; p.ub(2) = 0;
        else
            p.lb(3) = 1;
        end
        q = propagate_binary_product(p);
        testCase.verifyEqual(q.feasible,1);
        if fixed == 1
            testCase.verifyEqual(q.lb(3),1);
        elseif fixed < 4
            testCase.verifyEqual(q.ub(3),0);
        else
            testCase.verifyEqual(q.lb(1:2),[1;1]);
        end
    end
end

function test_implied_integer_coefficients(testCase)
p.feasible = 1; p.K.f = 1;
p.integer_variables = 1; p.binary_variables = [];
rows = {[0 1 -2],[0 .5 -1],[.5 1 -1],[0 2 -1],[3 -2 1]};
expected = {[],[],[],2,2};
for i = 1:numel(rows)
    p.F_struc = sparse(rows{i});
    q = presolve_impliedintegers_from_equalities(p);
    testCase.verifyEqual(q.implied_integers(:),expected{i}(:));
end

function test_periodic_objective_dependence(testCase)
p.c = [0;1]; p.Q = sparse(2,2);
p.lb = [-10;-Inf]; p.ub = [10;Inf];
p.monomtable = speye(2); p.F_struc = sparse(0,3);
p.evalMap = {struct('variableIndex',1,'properties',struct('periodic',2*pi))};
q = presolve_boundargument_periodic(p);
testCase.verifyEqual(q.ub(1),-10+2*pi,'AbsTol',1e-12);
p.c(1) = -1;
q = presolve_boundargument_periodic(p);
testCase.verifyEqual(q.ub,p.ub);
p.c(1) = 0; p.Q(1,1) = 1;
q = presolve_boundargument_periodic(p);
testCase.verifyEqual(q.ub,p.ub);
p.Q = sparse([0 1;1 0]);
q = presolve_boundargument_periodic(p);
testCase.verifyEqual(q.ub,p.ub);

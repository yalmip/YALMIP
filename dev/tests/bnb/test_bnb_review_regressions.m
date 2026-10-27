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

function test_fractional_knapsack_rhs(testCase)
yalmip('clear');
x = binvar(3,1);
ops = sdpsettings('solver','bnb','bnb.solver','linprog', ...
                  'bnb.maxiter',20,'verbose',0);
sol = optimize(0.6*sum(x) <= 1.5,-sum(x),ops);
testCase.verifyEqual(sol.problem,0);
testCase.verifyEqual(value(sum(x)),2,'AbsTol',1e-7);

function test_objective_constant_in_projected_sdp_cut(testCase)
p.c = [0;1]; p.f = -10;
p.lb = [0;0]; p.ub = [1;2];
p.noninteger_variables = 2; p.UnitContinuousCost = true;
row = fix_continuous_in_row_at_best_possible([-0.5 1 1],p,-9);
% b=0, t=1 is feasible for both the SDP cut and the incumbent bound.
testCase.verifyGreaterThanOrEqual(row*[1;0;1],-1e-12);
testCase.verifyEqual(row,[0.5 1 0],'AbsTol',1e-12);

function test_objective_constant_in_sdp_bound(testCase)
p.c = [0;1]; p.f = -10;
p.lb = [0;0]; p.ub = [1;2];
p.binary_variables = 1;
p.binarycardinality.up = 1; p.binarycardinality.down = 0;
sdpCuts.F_struc = sparse([-0.5 1 1]);
% b=0, t=1 satisfies the cut and has objective -9.
testCase.verifyLessThanOrEqual(upper_from_sdpcuts(p,sdpCuts,-9),-9);

function test_objective_cover_preserves_feasible_points(testCase)
p.c = 5*ones(4,1); p.f = -10;
p.LinearBinaryPositiveCost = true;
p.knapsack.a = {}; p.knapsack.b = {};
p.knapsack.variables = {}; p.knapsack.type = [];
cuts = knapsack_add_cover_cut(p,ones(4,1),'gu',5);
for mask = 0:15
    x = bitget(mask,1:4)';
    if p.f + p.c'*x <= 5
        testCase.verifyGreaterThanOrEqual(full(cuts*[1;x]),-1e-12);
    end
end

function test_sdp_cuts_from_all_cones_and_binary_only_covers(testCase)
p.K.s = [1 1];
p.semidefinite{1}.F_struc = sparse([-2 1 1 0]);
p.semidefinite{2}.F_struc = sparse([-0.5 0 0 1]);
p.options = sdpsettings;
p.lb = [0;0;0]; p.ub = [2;1;1];
p.isbinary = [false;true;true];
p.noninteger_variables = [];
p.UnitContinuousCost = false;
[covers,binaryCuts,rawCuts] = create_sdpknapsack_cuts(p,[0;0;0],p,inf);
testCase.verifyEqual(rawCuts.K.l,2);
testCase.verifyEqual(binaryCuts.K.l,1);
testCase.verifyEqual(full(binaryCuts.F_struc(1,2)),0);
% z=2, b=0, c=1 satisfies both SDP cuts; a 0-1 cover on z can exclude it.
point = [1;2;0;1];
testCase.verifyGreaterThanOrEqual(full(rawCuts.F_struc*point),-1e-12);
if covers.K.l > 0
    testCase.verifyGreaterThanOrEqual(full(covers.F_struc*point),-1e-12);
end

function test_cutsdp_cover_with_general_integer(testCase)
yalmip('clear');
z = intvar(1); b = binvar(1);
F = [0 <= z <= 2, [z+b-1 1;1 1] >= 0];
ops = sdpsettings('solver','cutsdp','cutsdp.solver','intlinprog', ...
                  'cutsdp.maxiter',20,'verbose',0);
sol = optimize(F,b+0.01*z,ops);
testCase.verifyEqual(sol.problem,0);
testCase.verifyEqual([value(z) value(b)],[2 0],'AbsTol',1e-7);

function tests = test_bnb_onehot_levels
tests = functiontests(localfunctions);

function test_negated_onehot_level_recovery(testCase)
% Two independent three-level groups enter the SDP only by their levels.
yalmip('clear');
x = intvar(6,1);
first = x(1:3);
second = x(4:6);
level1 = -[1 2 3]*first;
level2 = -[1 2 3]*second;
matrix = [level1-1.5 0.25; 0.25 level2-0.5];
F = [-1<=x<=0,sum(first)>=-1,sum(second)>=-1,matrix>=0];
h = level1+level2;
ops = sdpsettings('solver','bnb','verbose',0,'savesolveroutput',1);
sol = optimize(F,h,ops);
testCase.verifyEqual(sol.problem,0);
testCase.verifyEqual(value(h),3,'AbsTol',1e-6);
testCase.verifyEqual(value(x),[0;-1;0;-1;0;0],'AbsTol',1e-6);
testCase.verifyGreaterThanOrEqual(min(eig(value(matrix))),-1e-7);

function test_threshold_coordinate_change_preserves_all_rows(testCase)
% Check every original row and objective over all 16 feasible choices.
p.c = [1;2;3;1;2;3];
p.Q = sparse(6,6);
p.K.f = 0; p.K.l = 2; p.K.s = 2;
p.lb = zeros(6,1); p.ub = ones(6,1);
p.binary_variables = 1:6;
p.variabletype = zeros(6,1);
p.evalMap = {};
p.nonlinear = [];
p.rowtype = zeros(1,2);
p.x0 = [];
linear = [1 -1 -1 -1 0 0 0;1 0 0 0 -1 -1 -1];
S0 = [2 .1;.1 2];
cone = zeros(4,6);
cone(1,1:3) = 1:3;
cone(4,4:6) = 1:3;
p.F_struc = sparse([linear;[S0(:) cone]]);
[q,groups] = reencodeOneHotLevelGroups(p,1:6);
testCase.verifyEqual(numel(groups),2);
testCase.verifyEqual(q.K.l,6);
p.x0 = [0;1;0;0;0;1];
[withInitialGuess,~] = reencodeOneHotLevelGroups(p,1:6);
testCase.verifyEqual(withInitialGuess.x0,[1;1;0;1;1;1]);
oldRows = [1:p.K.l,q.K.l+1:size(q.F_struc,1)];
for level1 = 0:3
    for level2 = 0:3
        x = zeros(6,1);
        y = zeros(6,1);
        if level1>0
            x(level1) = 1;
            y(1:level1) = 1;
        end
        if level2>0
            x(3+level2) = 1;
            y(4:3+level2) = 1;
        end
        testCase.verifyEqual(full(q.F_struc(oldRows,:)*[1;y]), ...
            full(p.F_struc*[1;x]),'AbsTol',1e-12);
        testCase.verifyEqual(q.c'*y,p.c'*x,'AbsTol',1e-12);
        testCase.verifyGreaterThanOrEqual( ...
            min(full(q.F_struc(p.K.l+1:q.K.l,:)*[1;y])),-1e-12);
    end
end

function tests = test_bnb_objective_counts
tests = functiontests(localfunctions);
end

function test_unique_cost_tuple_fixes_variables(testCase)
p=model([sqrt(2);sqrt(2);pi;pi],zeros(4,1),ones(4,1),0);
s=bnbInitObjectiveCounts(p,2*sqrt(2)+2*pi);
[q,~,info]=bnbPropagateObjectiveCounts(s,p,2.8,2.9);
testCase.verifyTrue(logical(q.feasible));
testCase.verifyEqual(info.fixed,4);
testCase.verifyEqual(q.lb,[1;1;0;0]);
testCase.verifyEqual(q.ub,q.lb);
end

function test_count_cut_without_individual_fixing(testCase)
p=model(sqrt(2)*ones(3,1),zeros(3,1),ones(3,1),0);
p.x0=0.5*ones(3,1);
s=bnbInitObjectiveCounts(p,3*sqrt(2));
[q,cuts,info]=bnbPropagateObjectiveCounts(s,p,1.4,1.5);
testCase.verifyEqual(info.fixed,0);
testCase.verifyEqual(q.lb,p.lb); testCase.verifyEqual(q.ub,p.ub);
testCase.assertNotEmpty(cuts);
testCase.verifyLessThan(full(cuts*[1;p.x0]),-0.1);
points=eye(3);
testCase.verifyGreaterThanOrEqual(full([ones(3,1) points]*cuts'),-1e-12);
end

function test_empty_window_and_cache_scope(testCase)
p=model(sqrt(2)*ones(3,1),zeros(3,1),ones(3,1),0);
s=bnbInitObjectiveCounts(p,3*sqrt(2));
q=bnbPropagateObjectiveCounts(s,p,1.5,1.8);
testCase.verifyFalse(logical(q.feasible));
wide=p; wide.ub(1)=2;
[q,cuts]=bnbPropagateObjectiveCounts(s,wide,1.4,1.5);
testCase.verifyEqual(q,wide); testCase.verifyEmpty(cuts);
[q,cuts]=bnbPropagateObjectiveCounts(s,p,1.4,4*sqrt(2));
testCase.verifyEqual(q,p); testCase.verifyEmpty(cuts);
end

function test_signed_integer_oracle_preserves_every_feasible_point(testCase)
rng(611,'twister'); n=6;
digits=dec2base(0:3^n-1,3,n)-'0';
for trial=1:30
    atoms=[0.5;sqrt(2);sqrt(3)];
    c=atoms(randi(3,n,1)).*sign(rand(n,1)-0.5)+1e-12*randn(n,1);
    lo=randi([-3 2],n,1); p=model(c,lo,lo+2,randn);
    points=digits+lo'; values=p.f+points*c; ordered=sort(values);
    cap=ordered(randi([500 numel(ordered)]));
    s=bnbInitObjectiveCounts(p,cap); testCase.assertNotEmpty(s);
    for k=1:15
        node=p; witness=points(randi(size(points,1)),:)'; fixed=rand(n,1)>0.5;
        node.lb(fixed)=witness(fixed); node.ub(fixed)=witness(fixed);
        node.x0=node.lb+(node.ub-node.lb).*rand(n,1);
        upper=ordered(1)+(cap-ordered(1))*rand;
        lower=ordered(1)+(upper-ordered(1))*rand;
        feasible=all(points>=node.lb' & points<=node.ub',2) & values>=lower & values<=upper;
        [q,cuts]=bnbPropagateObjectiveCounts(s,node,lower,upper);
        if ~q.feasible, testCase.verifyFalse(any(feasible)); end
        remaining=all(points>=q.lb' & points<=q.ub',2);
        testCase.verifyTrue(all(remaining(feasible)));
        if ~isempty(cuts)
            slack=[ones(nnz(feasible),1) points(feasible,:)]*cuts';
            testCase.verifyTrue(all(slack(:)>=-1e-10));
        end
    end
end
end

function p=model(c,lo,hi,f)
n=numel(c);
p=struct('c',c,'f',f,'Q',sparse(n,n),'evalMap',[], ...
    'variabletype',zeros(n,1),'K',struct('m',0),'IntegerCost',false, ...
    'integral_variables',1:n,'lb',lo,'ub',hi,'feasible',1,'x0',(lo+hi)/2);
end

function tests = test_operator_review_regressions
tests = functiontests(localfunctions);

function test_logsumexp(testCase)
x = sdpvar(2,1);
y = logsumexp(x);
[~,op] = logsumexp('callback',[],x);
for offset = [-1000 0 1000]
    z = offset+[0;log(2)];
    assign(x,z);
    testCase.verifyEqual(logsumexp(z),offset+log(3),'AbsTol',1e-12);
    testCase.verifyEqual(value(y),offset+log(3),'AbsTol',1e-12);
    testCase.verifyEqual(op.derivative(z),[1;2]/3,'AbsTol',1e-12);
    [lo,hi] = op.bounds(z,z+1);
    testCase.verifyEqual([lo hi],offset+log(3)+[0 1],'AbsTol',1e-12);
end
testCase.verifyEqual(logsumexp([-Inf;-Inf]),-Inf);
testCase.verifyEqual(logsumexp([Inf;0]),Inf);
testCase.verifyEqual(logsumexp([]),-Inf);

function test_blackbox_metadata(testCase)
x = sdpvar(1);
hull = @(L,U)deal([],[],[],[]);
for direction = {'increasing','decreasing'}
    y = blackbox(x,@(z)z,direction{1},'convexhull',hull);
    e = yalmip('extstruct',getvariables(y));
    op = e.arg{3};
    testCase.verifyEqual(op.monotonicity,direction{1});
    testCase.verifyEqual(op.convexhull,hull);
end

function test_interpolation_nonuniform(testCase)
x = sdpvar(1); t = sdpvar(1);
xi = [0;1;10]; yi = [0;2;5];
for direction = [-1 1]
    [F,op] = interp1_internal('graph',t,x,xi,direction*yi,'graph',[]);
    if direction == 1
        testCase.verifyEqual(op.convexity,'concave');
    else
        testCase.verifyEqual(op.convexity,'convex');
    end
    for i = 1:3
        assign(x,xi(i)); assign(t,direction*yi(i));
        testCase.verifyGreaterThanOrEqual(min(check(F)),-1e-12);
    end
    [~,reversed] = interp1_internal('graph',t,x,flipud(xi),direction*flipud(yi),'graph',[]);
    testCase.verifyEqual(reversed.convexity,op.convexity);
end

function test_interpolation_monotonicity(testCase)
x = sdpvar(1); t = sdpvar(1);
data = {[2;1;0],[0;1;2],[0;2;1]};
expected = {'decreasing','increasing','none'};
for i = 1:3
    [~,op] = interp1_internal('graph',t,x,[0;1;2],data{i},'graph',[]);
    testCase.verifyEqual(op.monotonicity,expected{i});
end

function test_sine_bounds(testCase)
checkPeriodicBounds(testCase,@sin,[-3 2.8],[-1 1]);

function test_cosine_bounds(testCase)
checkPeriodicBounds(testCase,@cos,[-5.5 .5],[-1 1]);

function checkPeriodicBounds(testCase,f,interval,expected)
x = sdpvar(1); y = f(x);
e = yalmip('extstruct',getvariables(y));
[~,op] = feval(e.fcn,'graph',e.var,e.arg{1:end-1});
for shift = (-4:4)*2*pi
    [lo,hi] = op.bounds(interval(1)+shift,interval(2)+shift);
    testCase.verifyEqual([lo hi],expected,'AbsTol',1e-12);
end
for left = linspace(-12,12,17)
    for width = [.1 1 3 6 8]
        [lo,hi] = op.bounds(left,left+width);
        values = f(linspace(left,left+width,101));
        testCase.verifyLessThanOrEqual(lo,min(values)+1e-12);
        testCase.verifyGreaterThanOrEqual(hi,max(values)-1e-12);
    end
end

function test_tangent_bounds(testCase)
x = sdpvar(1); y = tan(x);
e = yalmip('extstruct',getvariables(y));
[~,op] = feval(e.fcn,'graph',e.var,e.arg{1:end-1});
for shift = (-4:4)*pi
    [lo,hi] = op.bounds(-2+shift,-1+shift);
    testCase.verifyEqual([lo hi],[-Inf Inf]);
    [lo,hi] = op.bounds(-.2+shift,.3+shift);
    testCase.verifyEqual([lo hi],tan([-.2 .3]),'AbsTol',1e-12);
end
for limits = {[-pi/2 -1],[0 pi/2]}
    z = limits{1};
    [lo,hi] = op.bounds(z(1),z(2));
    testCase.verifyEqual([lo hi],[-Inf Inf]);
end

function test_sum_square_bounds(testCase)
[~,op] = sum_square('callback',[],sdpvar(2,1));
[lo,hi] = op.bounds([-3;-4],[-2;-1]);
testCase.verifyEqual([lo hi],[5 25]);
[lo,hi] = op.bounds([-3;1],[2;4]);
testCase.verifyEqual([lo hi],[1 25]);

function test_xexp_convexity(testCase)
[~,op] = xexp('callback',[],sdpvar(1));
testCase.assertClass(op.convexity,'function_handle');
testCase.verifyEqual(op.convexity(-4,-2),'concave');
testCase.verifyEqual(op.convexity(-2,3),'convex');
testCase.verifyEqual(op.convexity(-4,3),'none');

function test_pexp_interior_minimum(testCase)
[~,op] = pexp('callback',[],sdpvar(2,1));
[lo,hi] = op.bounds([.5;1],[2;1]);
testCase.verifyEqual(lo,exp(1),'AbsTol',1e-12);
testCase.verifyEqual(hi,.5*exp(2),'AbsTol',1e-12);
for y = [-2 0 1 3]
    [lo,hi] = op.bounds([.5;y],[2;y+1]);
    values = arrayfun(@(t)pexp([t;y]),linspace(.5,2,21));
    testCase.verifyLessThanOrEqual(lo,min(values)+1e-12);
    testCase.verifyGreaterThanOrEqual(hi,max(values)-1e-12);
end
[lo,hi] = op.bounds([0;1],[2;1]);
testCase.verifyEqual([lo hi],[0 Inf]);

function test_plog_unbounded_box(testCase)
[~,op] = plog('callback',[],sdpvar(2,1));
[lo,hi] = op.bounds([1;1],[1;Inf]);
testCase.verifyEqual(hi,Inf);
testCase.verifyLessThanOrEqual(lo,0);
[lo,hi] = op.bounds([1;1],[Inf;2]);
testCase.verifyEqual(lo,-Inf);
testCase.verifyGreaterThanOrEqual(hi,log(2));

function test_slogfrac_denominator_pole(testCase)
[~,op] = slogfrac('callback',[],sdpvar(2,1));
for limits = {[-1 1],[-1 0],[-sqrt(eps) 1]}
    z = limits{1};
    [lo,hi] = op.bounds([1;z(1)],[1;z(2)]);
    testCase.verifyEqual([lo hi],[-Inf Inf]);
end

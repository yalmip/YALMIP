function tests = test_sos_1
tests = functiontests({@test1});

function test1(testCase)
yalmip('clear');

ops{1} = sdpsettings('sos.cong',0,'sos.model',1,'verbose',0);
ops{2} = sdpsettings('sos.cong',1,'sos.model',2,'verbose',0);
ops{3} = sdpsettings('sos.cong',0,'sos.newton',0,'verbose',0,'sos.extlp',0);

sdpvar x s t;
F = (sos(1+x+(1-s)*x^2-s))+(sos(2+2*x+x^4-8*t))+(s>=0.49)+(s+t<=0.55);
obj= -s-t;

for i = 1:length(ops)
    fail = regresstest(F,obj,ops{i});
    testCase.assertTrue(fail == 0);
end

function fail  = regresstest(F,obj,ops,pv);

if nargin==3
    pv = [];
end

ops.sos.model = 1;
solvesos(F,obj,ops,pv);
obj1 = value(obj);
p1s = check(F(find(is(F,'sos'))));
p1e = check(F(find(~is(F,'sos'))));

ops.sos.model = 2;
solvesos(F,obj,ops,pv);
obj2 = value(obj);
p2s = check(F(find(is(F,'sos'))));
p2e = check(F(find(~is(F,'sos'))));

fail = 0;

if abs(obj1-obj2) > 1e-4
    fail = 1;
end
if any(p1s>1e-4)
   fail = 2;
end
if any(p2s>1e-4)
   fail = 2;
end
if any(p1e<-1e-4)
   fail = 2;
end
if any(p2e<-1e-4)
   fail = 2;
end
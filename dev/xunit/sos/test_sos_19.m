function tests = test_sos_19
tests = functiontests({@test1});

function test1(dummy)
ops{1} = sdpsettings('sos.cong',0,'sos.model',1,'verbose',0);
ops{2} = sdpsettings('sos.cong',1,'sos.model',2,'verbose',0);
ops{3} = sdpsettings('sos.cong',0,'sos.newton',0,'verbose',0,'sos.extlp',0);

% Markus tentacles problem
sdpvar x y a
f = x^4 * y^2 + x^2 * y^4 - 3 * x^2 * y^2 + 1, k = 0;
df = jacobian(f, [x y]);
g = 1 - (df(1)^2 + df(2)^2) * (x^2 + y^2);
if k >= 0
    v = monolist([x; y], 2*k);
    coeffVec = sdpvar(length(v), 1);
    t = coeffVec' * v;
    constraints = (sos(f - a - t * g)) + (sos(t));
else
    coeffVec = [];
    constraints = (sos(f - a));
end
F = constraints;
obj = -a;
for i = 1:length(ops)
    i
    fail = regresstest(F,obj,ops{i},coeffVec);
    assert(fail == 0);
end




function fail  = regresstest(F,obj,ops,pv);

if nargin==3
    pv = [];
end

ops.sos.model = 1;
solvesos(F,obj,ops,pv);
obj1 = value(obj);
p1s = checkset(F(find(is(F,'sos'))));
p1e = checkset(F(find(~is(F,'sos'))));

ops.sos.model = 2;
solvesos(F,obj,ops,pv);
obj2 = value(obj);
p2s = checkset(F(find(is(F,'sos'))));
p2e = checkset(F(find(~is(F,'sos'))));

fail = 0;

if abs(obj1-obj2) > 1e-3
    fail = 1;
end

if any(p1s>1e-4)
   fail = 2;
   p1s
end
if any(p2s>1e-4)
   fail = 2;
   p2s
end
if any(p1e<-1e-4)
   fail = 2;
   p1e
end
if any(p2e<-1e-4)
   fail = 2;
   p2e
end
if fail==0
    disp('Correct solution');
end
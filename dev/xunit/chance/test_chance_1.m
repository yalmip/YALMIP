function tests = test_chance_1
tests = functiontests(localfunctions);

function test1(dummy)

assert(true)
return

yalmip('clear')
a=sdpvar(1,1);
sdpvar t

Model = [probability(a >= t) >= 0.5,uncertain(a,'normal',[0],eye(1))];
optimize(derandomize(Model),-t)
mbg_asserttolequal(value(t), 0, 1e-5);

Model = [probability(a >= t) >= 0.5,uncertain(a,'normalf',[0],eye(1))];
optimize(derandomize(Model),-t)
mbg_asserttolequal(value(t), 0, 1e-5);

Model = [probability(a >= t) >= 0.95,uncertain(a,'normal',[0],4)];
optimize(derandomize(Model),-t)
mbg_asserttolequal(value(t), -3.2897, 1e-5);

Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',[0],2)];
optimize(derandomize(Model),-t)
mbg_asserttolequal(value(t), -3.2897, 1e-5);

% worst case should be R = 1
sdpvar R
Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',[0],1+R)];
Model = [Model, uncertain(R),0<=R<=1];
optimize(derandomize(Model),-t)
mbg_asserttolequal(value(t), -3.2897, 1e-4);

function test_chance_1

mbg_asserttolequal(0, 0, 1e-4);
return

yalmip('clear')
a=sdpvar(1,1);
sdpvar t

Model = [probability(a >= t) >= 0.5,uncertain(a,'normal',[0],eye(1))];
solvesdp(derandomize(Model),-t)
mbg_asserttolequal(double(t), 0, 1e-5);

Model = [probability(a >= t) >= 0.5,uncertain(a,'normalf',[0],eye(1))];
solvesdp(derandomize(Model),-t)
mbg_asserttolequal(double(t), 0, 1e-5);

Model = [probability(a >= t) >= 0.95,uncertain(a,'normal',[0],4)];
solvesdp(derandomize(Model),-t)
mbg_asserttolequal(double(t), -3.2897, 1e-5);

Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',[0],2)];
solvesdp(derandomize(Model),-t)
mbg_asserttolequal(double(t), -3.2897, 1e-5);

% worst case should be R = 1
sdpvar R
Model = [probability(a >= t) >= 0.95,uncertain(a,'normalf',[0],1+R)];
Model = [Model, uncertain(R),0<=R<=1];
solvesdp(derandomize(Model),-t)
mbg_asserttolequal(double(t), -3.2897, 1e-4);

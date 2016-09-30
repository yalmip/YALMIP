function test_chance_2

yalmip('clear')
sdpvar a b t
sdpvar t

Model = [probability(a >= t) >= 0.5,
         probability(b >= t) >= 0.5, 
         uncertain(a,'normal',[0],eye(1)),
         uncertain(b,'normal',[0],eye(1))];
solvesdp(derandomize(Model),-t)
mbg_asserttolequal(double(t), 0, 1e-5);


yalmip('clear')
sdpvar a b t m
sdpvar t

Model = [probability(a >= t) >= 0.5,
         probability(b >= t) >= 0.5, 
         uncertain(a,'normal',[0],eye(1)),
         uncertain(b,'normal',m,eye(1)),
         uncertain(m,'normal',0,eye(1)), ];
P=optimizer(derandomize(Model),-t,sdpsettings('solver','mosek'),m,t)
Q=sample(P,20)

sdpvar w wmean
optimize(derandomize([uncertain(w,'normal',wmean,1), probability(w >= 0) >= 0.9]),wmean,sdpsettings('solver','cplex','debug',1))
mbg_asserttolequal(double(wmean), 1.2816, 1e-3);

clear all;yalmip('clear');
sdpvar w(2,1) wmean(2,1) gamma
uu = [];
for gamma = .1:.05:.9
    Model = derandomize([uncertain(w,'normal',wmean,1), probability(sum(w) >= 0) >= gamma]);
    optimize(Model,wmean'*wmean-gamma,sdpsettings('solver','','debug',1));
    uu = [uu value(wmean'*wmean-gamma)];
end
end



mbg_asserttolequal(double(wmean), 1.2816, 1e-3);

nnz([1 1]*[repmat(value(wmean),1,1000) + randn(2,1000)]>=0)

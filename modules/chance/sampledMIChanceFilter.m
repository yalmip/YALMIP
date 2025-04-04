function newConstraint =  sampledMIChanceFilter(b,c,distribution,gamma,w,options);
N = options.chance.N;
W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
d = binvar(N,1);
newConstraint = [implies(d,-b-c'*W <= 0), sum(d) >= N*(1-gamma)];
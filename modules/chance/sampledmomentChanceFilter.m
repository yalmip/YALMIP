function newConstraint =  sampledmomentChanceFilter(b,c,distribution,confidencelevel,w,options);
N = options,chance.N;
W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
d.parameters{2} = mean(W,2);
d.parameters{3} = cov(W');;
newConstraint = momentChanceFilter(b,c,d,confidencelevel,w,options);
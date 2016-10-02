function newConstraint =  sampledmomentChanceFilter(b,c,distribution,confidencelevel,w);
W = [];for i = 1:1000;W = [W dataSampler(distribution,size(w))];end
d.parameters{2} = mean(W,2);
d.parameters{3} = cov(W');;
newConstraint = momentChanceFilter(b,c,d,confidencelevel);
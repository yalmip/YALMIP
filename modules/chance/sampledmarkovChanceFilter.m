function newConstraint =  sampledmarkovChanceFilter(b,c,distribution,confidencelevel,w);
W = [];for i = 1:100;W = [W dataSampler(distribution,size(w))];end
alpha = sdpvar(1);
s = sdpvar(1,100);
newConstraint = [0 <= alpha, -b-c'*W+alpha <= s, 0 <= s, sum(s)/length(s) <= alpha*(1-confidencelevel)]

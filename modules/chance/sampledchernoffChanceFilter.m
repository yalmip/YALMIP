function newConstraint =  sampledchernoffChanceFilter(b,c,distribution,confidencelevel,w);
N = 250;
W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
alpha = sdpvar(1);
e = pexpsum([repmat(alpha,1,N);-b-c'*W])/N;
newConstraint = [0 <= alpha, e <= alpha*(1-confidencelevel)];
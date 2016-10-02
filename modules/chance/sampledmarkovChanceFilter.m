function newConstraint =  sampledmarkovChanceFilter(b,c,distribution,confidencelevel,w,options);
N = options.chance.N;
if strcmpi(func2str(distribution.name),'random') && strcmpi(distribution.parameters{1},'data')
    W = distribution.parameters{2};
    s = sdpvar(1,size(W,2));
else
    W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
    s = sdpvar(1,N);
end
alpha = sdpvar(1);
newConstraint = [0 <= alpha, -b-c'*W+alpha <= s, 0 <= s, sum(s)/length(s) <= alpha*(1-confidencelevel)]

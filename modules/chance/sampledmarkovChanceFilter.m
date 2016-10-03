function newConstraint =  sampledmarkovChanceFilter(b,c,distribution,confidencelevel,w,options);
if strcmpi(func2str(distribution.name),'random') && strcmpi(distribution.parameters{1},'data')
    W = distribution.parameters{2};
    N = size(W,2);
    s = sdpvar(1,N);
else
    N = options.chance.N;
    W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
    s = sdpvar(1,N);
end
alpha = sdpvar(1);
newConstraint = [0 <= alpha, -b-c'*W+alpha <= s, 0 <= s, sum(s)/N <= alpha*(1-confidencelevel)];

function newConstraint =  sampledmomentChanceFilter(b,c,distribution,confidencelevel,w,options);
if strcmpi(func2str(distribution.name),'random') && strcmpi(distribution.parameters{1},'data')
    W = distribution.parameters{2};
    N = size(W,2);   
else
    N = options.chance.N;
    W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
    s = sdpvar(1,N);
end
d.parameters{2} = mean(W,2);
d.parameters{3} = cov(W');;
newConstraint = momentChanceFilter(b,c,d,confidencelevel,w,options);
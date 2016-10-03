function newConstraint =  sampledchebyshevChanceFilter(b,c,distribution,confidencelevel,w,options);
N = options.chance.N;
W = [];for i = 1:N;W = [W dataSampler(distribution,size(w))];end
alpha = sdpvar(1);
s = sdpvar(1,N);
if isa(confidencelevel,'sdpvar')
    sdpvar t
    newConstraint = [-(b + c'*W) + alpha <= s, 0 <= s, norm_callback(s) <= alpha*sqrt(N)*sqrtm(1-confidencelevel)];
else
    newConstraint = [-(b + c'*W) + alpha <= s, 0 <= s, norm(s) <= alpha*sqrt(N)*sqrtm(1-confidencelevel)];
end

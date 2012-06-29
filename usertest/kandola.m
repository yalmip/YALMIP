tic
n=15;
K = randn(n);K = K*K';
t = sdpvar(1,1);
P = sdpvar(n,n);
bound = sdpvar(n*n,1);
s = sdpvar(1,1)
lambda = 2;
F = lmi(P>0)+lmi(cone(vec(P-K),t))+lmi(bound>P(:))+lmi(P(:)>-bound)+lmi([s t;t 1]>0)
toc
solvesdp(F,[], s +lambda*sum(bound))
P = double(P); P(1e-4>abs(P))=0
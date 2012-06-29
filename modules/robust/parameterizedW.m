function W = parameterizedW(G,H,x)

ng = size(G,1);
mg = size(H,2);

K = full(sparse(1:ng*mg,reshape(reshape(1:ng*mg,ng,mg)',ng*mg,1),1));
x = recover(depends(H));
n = length(x);

Hbase = getbase(H);
for i = 1:n
    Hi{i} = full(getbasematrix(H,getvariables(x(i))));
end

for i = 1:n
    S{i} = kron(eye(ng),Hi{i})*K + kron(Hi{i},eye(ng));
end

SS = [blkdiag(S{1},S{2},S{3});S{2} S{1} 0*S{3};S{3} 0*S{2} S{1};0*S{1} S{3} S{2}];
Snull = null(SS);
Snull(abs(Snull)<1e-6) = 0;
w = sdpvar(size(Snull,2),1);
Sw = Snull*w;


W = sdpvar(size(G,1),size(H,2));
top = 1;
for i = 1:n
    W = W + reshape(Sw(top:top+ng*mg-1),ng,mg)*x(i);
    top = top + ng*mg;
end

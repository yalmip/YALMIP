function cost = computecost(f,c,Q,x,p);

cost = f+c'*x+x'*Q*x;
if ~isequal(p.K.m,0)
    X = p.F_struc(end-p.K.m^2+1:end,:)*[1;x];
    X = reshape(X,p.K.m,p.K.m);
    cost = cost - sum(log(real(eig(X))));
end
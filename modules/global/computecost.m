function cost = computecost(f,c,Q,x,p);

cost = f+c'*x+x'*Q*x;
if ~isequal(p.K.m,0)
    top = size(p.F_struc,1)-sum(p.K.m.^2)+1;
    for i = 1:length(p.K.m)
        X = p.F_struc(top:top + p.K.m(i)^2-1,:)*[1;x];
        X = reshape(X,p.K.m(i),p.K.m(i));
        cost = cost - sum(real(log(real(eig(X)))));
        top = top + p.K.m(i)^2;
    end
end
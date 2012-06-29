function f = invsathub(p,lambda)

absp = abs(p);
f = lambda*absp;
case2=find(absp>lambda);
f(case2) = -0.25*(p(case2).^2-6*lambda*absp(case2)+lambda^2);
case3=find(absp>3*lambda);
f(case3) = 2*lambda^2;


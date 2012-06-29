function [y,F] = pwasos2(fi,xi,x)

fi = fi(:);
xi = xi(:);

n = length(fi)-1;
lambda = sdpvar(n+1,1);
reg    = binvar(n,1);

F = [sum(lambda)==1,lambda>0,sum(reg)==1];

F = [F,lambda(0+1) <= reg(1)];
for l =1:n-1
    F = [F,lambda(l+1)-reg(l)-reg(l+1) < 0];
end
F = [F,lambda(end)<reg(end)];
F = [F,lambda'*xi(:) == x];
y = lambda'*fi;


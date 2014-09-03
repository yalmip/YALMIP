function [p,A,b] = pwamodel(f,x)

t = sdpvar(1);
F = ([f <= t]);
[A,b] = pwamodel(F,[x;t]);
B = A(:,end);
A = A(:,1:end-1);

% Ax + Bt <= b
for i = 1:length(B)
    if B(i)~=0
        A(i,:) = A(i,:)/abs(B(i));
        b(i) = b(i)/abs(B(i));
    else
    end
end
p = max(A*x-b);
b = -b;

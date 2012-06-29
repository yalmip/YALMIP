function A = adjoint(X)
% ADJOINT Computes adjoint matrix
%
% A = ADJOINT(X)
%
% Brute-force implementation

[n,m] = size(X);
if n~=m
    error('Matrix must be square');
end

A = [];
if n == 1
    A = X;
    return
end

% Ugly brute-force
for i = 1:n
    temp = [];
    noti = setdiff(1:n,i);
    for j = 1:n
        notj = setdiff(1:n,j);
        temp = [temp det(X(noti,notj))*((-1)^(i+j))];
    end
    A = [A;temp];
end
A = A';


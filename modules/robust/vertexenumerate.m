function V = vertexenumerate(A,b,x0)
A = full(A);
b = full(b);
% Shift origin
b = b - A*x0;
if any(b<0)
    error('Shift failure in vertexenumerate')
end
n = size(A,2);
m = size(A,1);
% Scale b
A = A ./ repmat(b,1,n);
% Lift and find rays
k = convhulln([A;zeros(1,n)]);
V = [];
% Find solutions
rhs = ones(n,1);
for j = 1:size(k,1)
    V = [V A(k(j,:),:)\rhs];
end
% Map back
V = V + repmat(x0,1,size(V,2));
V = unique(V','rows')';
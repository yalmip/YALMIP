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
if n == 1
    % Silly
    V = [1/max(A) 1/min(A)];
    return
end
% Lift and find rays
k = convhulln([A;zeros(1,n)]);
if any(k(:)>size(A,1))
    error('Polytope is unbounded');
end
V = [];
% Find solutions
rhs = ones(n,1);
for j = 1:size(k,1)
    V = [V A(k(j,:),:)\rhs];
end
% Map back
V = V + repmat(x0,1,size(V,2));
V = unique(V','rows')';
function tbasis = createImagBasis(n)
nvari = n*(n+1)/2+(n*(n+1)/2-n);
tbasis = spalloc(n^2,1+nvari,2);
l = 2;
an_empty = spalloc(n,n,2);
Y = reshape(1:n^2,n,n);
Y = tril(Y);
Y = (Y+Y')-diag(sparse(diag(Y)));
[uu,oo,pp] = unique(Y(:));
BasisReal = sparse(1:n^2,pp+1,1);

% Remove diagonal, flip signs, and complexify
BasisImag = BasisReal;
BasisImag(:,1) = [];
index = [1];
for i = 1:n-1
    index = [index index(end)+n-i+1];
end
BasisImag(:,index) = [];
[i,j,s] = find(BasisImag);
s(1:2:end)=-1;
BasisImag = sparse(i,j,s*sqrt(-1),size(BasisImag,1),size(BasisImag,2));

tbasis = [BasisReal BasisImag];
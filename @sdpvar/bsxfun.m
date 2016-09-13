function out=bsxfun(f,x,y)
%BSXFUN

dim1 = size(x);
dim2 = size(y);

X = x;
Y = y;

for i = 1:length(dim1)
    if dim1(i) < dim2(i)
        d = ones(1,length(dim1));
        d(i) = dim2(i);
        X = repmat(X,d);
    elseif dim1(i) > dim2(i)
        d = ones(1,length(dim1));
        d(i) = dim1(i);
        Y = repmat(Y,d);
    end
end

% And compute
out = f(X,Y);